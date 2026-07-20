#!usr/bin/env bash

output=BATCH1
FQ_I1=SAMPLE_I1.fastq.gz
FQ_I2=SAMPLE_I2.fastq.gz
FQ_R1=SAMPLE_R1.fastq.gz
FQ_R2=SAMPLE_R2.fastq.gz

# generate reference barcode files
preprocess_barcodes.indel.pl -i barcodes -o $output -b I5 I7 T5 T7 -bf ./ -hd1 2 -hd2 2 -s spike-in

# demultiplexing
get_fixed.fq.indel.pl \
	-i $FQ_I1 $FQ_I2 $FQ_R1 $FQ_R2 \
	-b1 $output.demx.barcode.barcodes \
	-b2 $output.demx.barcode.samples \
	-s $output.demx.spikein.barcodes \
	-t 1 \
	-o $output
