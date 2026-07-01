#!/usr/bin/env perl
use Getopt::Long;
use Array::Utils qw(:all);

GetOptions(
	'input|i=s'			=> \$input,
	'qctable1|qc1=s'			=> \$qctable1,
	'minfrag|min=s'			=> \$minfrag,
	'maxfrag|max=s'			=> \$maxfrag,
	'nucleosome_signal|ns=s'			=> \$nucleosome_signal,
	'qctable2|qc2=s'			=> \$qctable2,
	'tss_enrich|te=s'			=> \$tss_enrich,
	'output|o=s'			=> \$output,
	'help|h!'			=> \$help,
);

if($help){
print "\nFiltering fragment: fragments count > 1000 & < 100000, TSS enrichment sore > 2, nucleosome signal < 4";
print "usage: $0\n";
print " -i: input fragment bed file, required\n";
print " -qc1: qc table 1 (output of 'fragment.qc.1.pl'), required\n";
print " -qc2: qc table 2 (output of 'ATACCellTSS'), required\n";
print " -o: output fragment file, default 'test'\n";
print " -h: help\n";

print "
 QC parameters:
  -min: minimal number of fragments in cells to keep, default '1000'
  -max: maximal number of fragments in cells to keep, default '100000'
  -ns: maximal nucleosome signal of cells to keep, default '4'
  -te: minimal TSS enrichment score of cells to keep, default '2'

";
print "Output format: cellbarcode\tALL_fragment_count.1-1000\tfragment_count.1-147\tfragment_count.147-294\tnucleosome_signal\tMT_ratio\tTSS_enrichment_score\n\n";

}elsif($input and $qctable1 and $qctable2){
	$output||="test";
	$minfrag||=1000;
	$maxfrag||=100000;
	$nucleosome_signal||=4;
	$tss_enrich||=2;
	
	## fragments count > 1000 & < 100000, nucleosome signal < 4
	open IN,"$qctable1" or die;
	my %qc1;
	my %qc3;
	while(<IN>){
		chomp;
		my @arr = split;
		$qc1{$arr[0]} = "$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]" if $arr[1] > $minfrag and $arr[1] < $maxfrag and $arr[4] < $nucleosome_signal;
		$qc3{$arr[0]} = "$arr[1]\t$arr[2]\t$arr[3]\t$arr[4]\t$arr[5]";
	}
	close IN;
	
	## TSS enrichment sore > 2
	open IN,"$qctable2" or die;
	my %qc2;
	while(<IN>){
		chomp;
		my @arr = split;
		$qc2{$arr[0]} = "$qc1{$arr[0]}\t$arr[1]" if $arr[1] > $tss_enrich and exists $qc1{$arr[0]};
		$qc3{$arr[0]} = "$qc3{$arr[0]}\t$arr[1]";
	}
	close IN;
	
	# out
	open IN,"$input" or die;
	open OUT,">$output" or die;
	while(<IN>){
		chomp;
		my @arr = split;
		print OUT "$_\n" if exists $qc2{$arr[3]};
	}
	close OUT;
	
	open OUT2,">$output.metadata" or die;
	foreach my $key (keys %qc3){
		print OUT2 "$key\t$qc3{$key}\n";
	}
	
}else{
	print "Didn't specify input files!\n";
	exec("$0 -h")
}


