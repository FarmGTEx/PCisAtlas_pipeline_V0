#!/usr/bin/env perl
use Getopt::Long;
use Array::Utils qw(:all);

GetOptions(
	'input|i=s'			=> \$input,
	'outfrag|o=s'			=> \$outfrag,
	'outstat|s=s'			=> \$outstat,
	'help|h!'			=> \$help,
);

if($help){
print "\nKeep fragments with length < 1000, and calculate the nucleosome signal";
print "usage: $0\n";
print " -i: input bed file (from 'bedtools intersect -c')\n";
print " -o: output fragment file, default 'test.frg'\n";
print " -s: output statistics of fragment file, default 'test.stat'\n";
print " -h: help\n";
print "Output format: cellbarcode\tALL_fragment_count.1-1000\tfragment_count.1-147\tfragment_count.147-294\tnucleosome_signal\tMT_ratio\n\n";

}elsif($input){
	$outfrag||="test.frg";
	$outstat||="test.stat";
	
	## filtering frags >= 1000 and generating cell hash
	if($input =~ /gz$/){
		open IN,"unpigz -dkc -p 1 $input|" or die;
	}else{
		open IN,"$input" or die;
	}
	open OUT,">$outfrag" or die;
	my %hash;
	my %mtfrags;
	while(<IN>){
		chomp;
		my @arr = split;
		if($arr[0] eq 'chrMT'){
			$mtfrags{$arr[3]}++;
			next;
		}
		my $fraglen = $arr[2] - $arr[1];
		if($fraglen < 147){
			$hash{$arr[3]}[0]++;
			$hash{$arr[3]}[1]++;
			print OUT "$_\n";
		}elsif($fraglen >= 147 and $fraglen < 294){
			$hash{$arr[3]}[0]++;
			$hash{$arr[3]}[2]++;
			print OUT "$_\n";
		}elsif($fraglen >= 294 and $fraglen < 1000){
			$hash{$arr[3]}[0]++;
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	
	open OUT,">$outstat" or die;
	foreach my $key (keys %hash){
		my $ratio = 999999;
		$ratio = $hash{$key}[2]/$hash{$key}[1] if $hash{$key}[1] != 0;
		my $mtratio = 0;
		$mtratio = $mtfrags{$key}/$hash{$key}[0] if exists $mtfrags{$key};
		print OUT "$key\t$hash{$key}[0]\t$hash{$key}[1]\t$hash{$key}[2]\t$ratio\t$mtratio\n";
	}
	close OUT;
	
}else{
	print "Didn't specify input files!\n";
	exec("$0 -h")
}


