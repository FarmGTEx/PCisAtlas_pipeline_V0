#!/usr/bin/env perl
use Getopt::Long;
use Array::Utils qw(:all);

GetOptions(
	'fragment|f=s'			=> \$fragment,
	'qctable|q=s'			=> \$qctable,
	'doublet|d=s'			=> \$doublet,
	'usedbarcode|u=s'		=> \$usedbarcode,
	'summarys|s=s'			=> \$summarys,
	'output|o=s'			=> \$output,
	'help|h!'			=> \$help,
);

if($help){
print "\nFiltering doublets with the result from 'AMULET'";
print "usage: $0\n";
print " -f: fragment file, required\n";
print " -q: qc table (meta data), required\n";
print " -d: barcodes of doublets from 'AMULET', required\n";
print " -u: barcodes used as input into 'AMULET' (barcodes filtered out the low quality barcodes), required\n";
print " -s: summary file (from 'QC.plot.r'), required\n";
print " -o: output fragment file, default 'test'\n";
print " -h: help\n";
print "Output format: cellbarcode\tALL_fragment_count.1-1000\tfragment_count.1-147\tfragment_count.147-294\tnucleosome_signal\tMT_ratio\tTSS_enrichment_score\tfragments_in_peaks\tRFiP\tDoublet\n\n";

}elsif($doublet and $fragment and $qctable and $summarys and $usedbarcode){
	$output||="test";
	
	## generate doublets hash
	open IN,"$doublet" or die;
	my %doublets;
	while(<IN>){
		chomp;
		$doublets{$_} = 1;
	}
	close IN;
	
	## generate barcodes hash
	open IN,"$usedbarcode" or die;
	my %usedbarcodes;
	while(<IN>){
		next if /^Cell/;
		chomp;
		my $key = (split)[0];
		$usedbarcodes{$key} = 1;
	}
	close IN;
	
	# output fragment and meta data
	open IN,"$fragment" or die;
	open OUT,">$output" or die;
	while(<IN>){
		chomp;
		my @arr = split;
		print OUT "$_\n" if not exists $doublets{$arr[3]};
	}
	close IN;
	close OUT;
	
	open IN,"$qctable" or die;
	open OUT,">$output.metadata" or die;
	while(<IN>){
		chomp;
		my @arr = split;
#		print OUT "$_\n" if not exists $doublets{$arr[0]};
		if(exists $doublets{$arr[0]}){
			print OUT "$_\tDoublet\n";
		}elsif(exists $usedbarcodes{$arr[0]}){
			print OUT "$_\tNonDoublet\n";
		}else{
			print OUT "$_\tNA\n";
		}
	}
	close IN;
	close OUT;
	
	# OUTPUT SUMMARY FILE
	my $double_num = readpipe("wc -l $doublet");
	$double_num = (split(" ",$double_num))[0];
	open IN,"$summarys" or die;
	open OUT,">$output.QC.summary" or die;
	while(<IN>){
		chomp;
		my @arr = split;
		print OUT "$_\n";
		if($arr[0] eq "Final_Cells"){
			my $final2 = $arr[1] - $double_num;
			print OUT "\nNumber of Multiplets\t$double_num\nFinal cells 2\t$final2\n";
		}
	}
	close IN;
	close OUT;
	
	
}else{
	print "Didn't specify input files!\n";
	exec("$0 -h")
}


