#!/usr/bin/env perl
use Getopt::Long;
use Array::Utils qw(:all);

GetOptions(
	'input|i=s'			=> \$input,
	'fragment|f=s'			=> \$fragment,
	'qctable|q=s'			=> \$qctable,
	'ratio_cut|r=s'			=> \$ratio_cut,
	'output|o=s'			=> \$output,
	'help|h!'			=> \$help,
);

if($help){
print "\nFiltering fragment: fragments count > 1000 & < 100000, TSS enrichment sore > 2, nucleosome signal < 4";
print "usage: $0\n";
print " -i: input bed file (output by 'bedtools intersect'), required\n";
print " -f: fragment file, required\n";
print " -q: qc table, required\n";
print " -r: threshold value of ratio of fragments in peaks, default '0.2'\n";
print " -o: output fragment file, default 'test'\n";
print " -h: help\n";
print "Output format: cellbarcode\tALL_fragment_count.1-1000\tfragment_count.1-147\tfragment_count.147-294\tnucleosome_signal\tMT_ratio\tTSS_enrichment_score\tfragments_in_peaks\tRFiP\n\n";

}elsif($input and $fragment and $qctable){
	$ratio_cut||=0.2;
	$output||="test";
	
	## generate and output file for QC
	open IN,"$input" or die;
	my %frg_all;
	my %frg_good;
	while(<IN>){
		chomp;
		my @arr = split;
		$frg_all{$arr[3]}++;
		$frg_good{$arr[3]}++ if $arr[4] > 0;
	}
	close IN;
	
	# output QC file (metadata)
	open IN,"$qctable" or die;
	open OUT,">$output.metadata" or die;
	my %ratios;
	while(<IN>){
		chomp;
		my @arr = split;
		$frg_good{$arr[0]} = 0 if not exists $frg_good{$arr[0]};
		my $ratio = $frg_good{$arr[0]}/$frg_all{$arr[0]};
		print OUT "$_\t$frg_good{$arr[0]}\t$ratio\n";
		$ratios{$arr[0]} = $ratio;
	}
	close IN;
	close OUT;
	
	# filter and output fragments
	open IN,"$fragment" or die;
	open OUT,">$output" or die;
	while(<IN>){
		chomp;
		my $barcode = (split)[3];
		print OUT "$_\n" if $ratios{$barcode} > $ratio_cut;
	}
	close IN;
	close OUT;
	
}else{
	print "Didn't specify input files!\n";
	exec("$0 -h")
}


