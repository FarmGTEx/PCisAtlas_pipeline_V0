#!/usr/bin/bash

# liftover
CRE_path=/path/to/CRE
chain_path=/path/to/liftOver_chain

liftOver ./02.Adipocyte.bed ./susScr11ToHg19.over.chain.gz 01.Adipocyte_pig_to_hg19.bed 02.Adipocyte_pig_to_hg19.unmapped.bed -minMatch=0.5 && \
awk -v FS='\t' -v OFS='\t' '{print $0,$1"-"$2"-"$3}' 01.Adipocyte_pig_to_hg19.bed > tmpfile_Adipocyte && \
mv tmpfile_Adipocyte 01.Adipocyte_pig_to_hg19.bed && \
liftOver 01.Adipocyte_pig_to_hg19.bed ./hg19ToSusScr11.over.chain.gz 03.Adipocyte_pig_to_hg19.re_pig.bed 04.Adipocyte_pig_to_hg19.re_pig.unmapped.bed -minMatch=0.5 && \
awk -v FS='\t' -v OFS='\t' '{print $0,$1"-"$2"-"$3}' 03.Adipocyte_pig_to_hg19.re_pig.bed > tmpfile_Adipocyte && \
mv tmpfile_Adipocyte 03.Adipocyte_pig_to_hg19.re_pig.bed && \
awk '{if($4==$6){print $5}}' 03.Adipocyte_pig_to_hg19.re_pig.bed | sed 's/-/\t/g' | bedtools sort -i - |awk -v FS='\t' -v OFS='\t' '{print $0,$1"-"$2"-"$3}' > 05.pig_Adipocyte.hg19_used.bed &


# run LDSC
ldsc.py \
  --h2-cts ./sumstats/PASS_HDL.sumstats.gz \
  --ref-ld-chr ./baseline_v1.2/baseline. \
  --out PASS_HDL \
  --ref-ld-chr-cts celltype.ldsc \
  --n-blocks 1000 \
  --w-ld-chr ./1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.
