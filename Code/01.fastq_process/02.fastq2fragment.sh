#!usr/bin/env bash
# processing fastq files and generate fragment and metadata for subsequent analysis

# 
NAME=SAMPLE
FQ_R1=BATCH1.demulti.SAMPLE_R1.fastq
FQ_R2=BATCH1.demulti.SAMPLE_R2.fastq
THREAD=20
ADAPTER_1=CTGTCTCTTA
ADAPTER_2=CTGTCTCTTA
AMULET=~/tools/AMULET-1.1
GENOME=~/refgenome/pig.fa
GENOMEINDEX=~/refgenome/pig.chromapindex
GTFFILE=~/refgenome/Sus_scrofa.Sscrofa11.1.101.chr.gtf
AUTOSOME=~/refgenome/pig_autosomes.txt
REPEAT=~/refgenome/pig.repeat.txt
BIN5K=~/refgenome/pig.5K
CHROMSIZE=~/refgenome/pig.sizes


##filtering fastq files by fastp
echo -e "Processing......filtering by fastp, data: `date`"
fastp -i $FQ_R1 -I $FQ_R2 \
 -o ${NAME}.filtered_R1.fastq.gz -O ${NAME}.filtered_R2.fastq.gz \
 -j ${NAME}.fastp.json -h ${NAME}.fastp.html \
 --thread $THREAD \
 -q 20 -u 20 -e 30 -l 20 \
 --adapter_sequence $ADAPTER_1 --adapter_sequence_r2 $ADAPTER_2 \
 --cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 \
 --cut_right --cut_right_window_size 3 --cut_right_mean_quality 20 \
 --dont_eval_duplication \
 --trim_poly_g


##align and filter
echo -e "Processing......chromap mapping and filtering, data: `date`"
chromap -t $THREAD -x $GENOMEINDEX -r $GENOME \
  -1 ${NAME}.filtered_R1.fastq.gz -2 ${NAME}.filtered_R2.fastq.gz \
  --min-read-length 20 \
  -l 2000 -q 30 --min-frag-length 30 \
  --SAM -o /dev/stdout \
  |sambamba view -q -t $THREAD -S -f bam /dev/stdin |sambamba sort -q -t $THREAD -n /dev/stdin -o ${NAME}.bam


#################### filtering procedure for fragments file
#remove duplicates
#fragments length < 1000, remove MT reads
#fragments count > 1000 & < 100000
#TSS enrichment sore > 2
#nucleosome signal < 4
#FRiP > 0.2
#remove doublets
##############


##generate fragment file
# removing duplicates and shifted the fragments by +4,-5
echo -e "Processing......generating fragment file, data: `date`"
bedtools bamtobed -bedpe -i ${NAME}.bam |awk 'BEGIN{FS="\t|:";OFS="\t"}{if($1==$4){print $1,$2+4,$6-5,$7}}' > ${NAME}.fragments.CB.shifted.bed
sort --parallel $THREAD -u ${NAME}.fragments.CB.shifted.bed |bedtools sort -i - > ${NAME}.fragments.CB.uniq.shifted.bed
QC.plot.fraglen.r -i ${NAME}.fragments.CB.uniq.shifted.bed -o ${NAME}.fragments.CB.uniq.shifted


## filer fragments, step 1
#fragments length < 1000, remove MT reads
echo -e "Processing......filer fragments (step 1), data: `date`"
fragment.qc.1.pl -i ${NAME}.fragments.CB.uniq.shifted.bed -o ${NAME}.fragments.CB.uniq.shifted.flt1.bed -s ${NAME}.fragments.CB.uniq.shifted.flt1.metadata


##TSS enrichment
# Computing overall TSS for the full dataset
# only test in snapatac2 version 2.3.0
echo -e "Processing......Computing TSS for the each barcode, data: `date`"
snapatac2.tsse.py --bed_file ${NAME}.fragments.CB.uniq.shifted.flt1.bed \
  --gtf_file $GTFFILE \
  --chromsizes_file $CHROMSIZE \
  --output_file ${NAME}.fragments.CB.uniq.shifted.flt1.tsse


## filer fragments (cell barcodes), step 2
#fragments count > 1000 & < 100000
#TSS enrichment sore > 2
#nucleosome signal < 4
echo -e "Processing......filter fragments (step 2), data: `date`"
fragment.qc.2.pl -i ${NAME}.fragments.CB.uniq.shifted.flt1.bed \
  -qc1 ${NAME}.fragments.CB.uniq.shifted.flt1.metadata \
  -qc2 ${NAME}.fragments.CB.uniq.shifted.flt1.tsse \
  -o ${NAME}.fragments.CB.uniq.shifted.flt2.bed \
  -min 1000 -max 100000 -ns 4 -te 2


## call peaks with MACS2 and filter cells by ratio of fragments in peaks > 0.2
echo -e "Processing......filter fragments by FRiP (step 3), data: `date`"
macs2 callpeak -f BEDPE -t ${NAME}.fragments.CB.uniq.shifted.flt2.bed -g 2367939388 -n ${NAME}.fragments.CB.uniq.shifted.flt2 \
  --shift 100 --extsize 200 --nomodel --call-summits --nolambda --keep-dup all --qval 1e-2
awk '{print $1,$2,$3,$1"_"$2"_"$3}' OFS="\t" ${NAME}.fragments.CB.uniq.shifted.flt2_peaks.narrowPeak > ${NAME}.fragments.CB.uniq.shifted.flt2_peaks.narrowPeak.bed
bedtools intersect -a ${NAME}.fragments.CB.uniq.shifted.flt1.bed -b ${NAME}.fragments.CB.uniq.shifted.flt2_peaks.narrowPeak.bed -c > ${NAME}.fragments.CB.uniq.shifted.flt1.CB_peak.bed
fragment.qc.3.pl -i ${NAME}.fragments.CB.uniq.shifted.flt1.CB_peak.bed -f ${NAME}.fragments.CB.uniq.shifted.flt2.bed -q ${NAME}.fragments.CB.uniq.shifted.flt2.bed.metadata -r 0.2 -o ${NAME}.fragments.CB.uniq.shifted.flt3.bed


## QC summary and plot
echo -e "Processing......QC plot 1, data: `date`"
QC.plot.r -i ${NAME}.fragments.CB.uniq.shifted.flt3.bed.metadata -f 1000 -r 0.2 -t 2 -n 4 -o ${NAME}.fragments.CB.uniq.shifted.flt3
nohup QC.plot.fraglen.r -i ${NAME}.fragments.CB.uniq.shifted.flt3.bed -o ${NAME}.fragments.CB.uniq.shifted.flt3 &> ${NAME}.nohup.qc1 &


## Doublets detection
echo -e "Processing......doublets detection with AMULET, data: `date`"
awk -v OFS='\t' '{print $0,1}' ${NAME}.fragments.CB.uniq.shifted.flt3.bed |bgzip > ${NAME}.fragments.CB.uniq.shifted.flt3.bed.tsv.gz && tabix -p bed ${NAME}.fragments.CB.uniq.shifted.flt3.bed.tsv.gz
awk -v OFS=',' '{if($2 > 1000 && $2 < 100000 && $5 < 4 && $7 > 2 && $9 > 0.2){print $1,1}}' ${NAME}.fragments.CB.uniq.shifted.flt3.bed.metadata |sed '1i\barcode,is__cell_barcode' > ${NAME}.fragments.CB.uniq.shifted.flt3.bed.metadata.csv
mkdir -p ${NAME}.doublet
AMULET.sh --maxinsertsize 1000 ${NAME}.fragments.CB.uniq.shifted.flt3.bed.tsv.gz ${NAME}.fragments.CB.uniq.shifted.flt3.bed.metadata.csv $AUTOSOME $REPEAT ${NAME}.doublet $AMULET
fragment.qc.4.pl -f ${NAME}.fragments.CB.uniq.shifted.flt3.bed -q ${NAME}.fragments.CB.uniq.shifted.flt3.bed.metadata -d ${NAME}.doublet/MultipletBarcodes_01.txt -u ${NAME}.doublet/OverlapSummary.txt -s ${NAME}.fragments.CB.uniq.shifted.flt3.QC.summary -o ${NAME}.fragments.CB.uniq.shifted.flt4.bed
Stat.fragments.r -i ${NAME}.fragments.CB.uniq.shifted.flt4.bed.metadata -s ${NAME}.fragments.CB.uniq.shifted.flt4.bed.QC.summary
sed -i 's/ tab /\t/' ${NAME}.fragments.CB.uniq.shifted.flt4.bed.QC.summary


## plot final
echo -e "Processing......QC plot (filnal), data: `date`"
nohup grep 'NonDoublet' ${NAME}.fragments.CB.uniq.shifted.flt4.bed.metadata > ${NAME}.fragments.CB.uniq.shifted.flt4.bed.metadata.tmp && \
  QC.plot.r -i ${NAME}.fragments.CB.uniq.shifted.flt4.bed.metadata.tmp -f 1000 -r 0.2 -t 2 -n 4 -o ${NAME}.fragments.CB.uniq.shifted.flt4 \
  &> ${NAME}.nohup.qc2 &
nohup QC.plot.fraglen.r -i ${NAME}.fragments.CB.uniq.shifted.flt4.bed -o ${NAME}.fragments.CB.uniq.shifted.flt4 &> ${NAME}.nohup.qc3 &

echo -e "All done!!! , data: `date`"
