#!/usr/bin/env python
# -*- coding:utf-8 -*-
#ref: https://kzhang.org/SnapATAC2/reference/_autosummary/snapatac2.pp.import_data.html
import argparse
import snapatac2
import scanpy as sc
import pandas as pd

def main(args):
	with open(args.chromsizes_file, 'r') as f:
		dict_chromsize = {}
		for line in f:
			key, value = line.strip().split()
			dict_chromsize[key] = int(value)
	
	data = snapatac2.pp.import_data(
		args.bed_file,
		gene_anno=args.gtf_file,
		chrom_size=dict_chromsize,
		sorted_by_barcode=False,
		min_num_fragments=0,
		min_tsse=0,
		low_memory=False
	)
	
	tsse_info = data.obs["tsse"]
	cell_names = data.obs_names
	data = pd.DataFrame({"cell_name": cell_names, "tsse_info": tsse_info})
	data.to_csv(args.output_file, sep="\t", index=False, header=False)


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Process some integers.')
	parser.add_argument('--bed_file', type=str, help='input fragment file')
	parser.add_argument('--gtf_file', type=str, help='path to annotation file (gtf/gff with "transcript")')
	parser.add_argument('--chromsizes_file', type=str, help='path to chromosome sizes file')
	parser.add_argument('--output_file', type=str, help='path to output file')

	args = parser.parse_args()
	main(args)
