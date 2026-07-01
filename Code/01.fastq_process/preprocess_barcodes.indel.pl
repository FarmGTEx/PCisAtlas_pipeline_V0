#!/usr/bin/env perl
use Getopt::Long;
use Array::Utils qw(:all);
use Math::Cartesian::Product;
use List::MoreUtils ':all';
use Text::Levenshtein qw(distance);

GetOptions(
	'input|i=s'				=> \$input,
	'output|o=s'			=> \$output,
	'barcode_folder|bf=s'	=> \$barcode_folder,
	'barcodes|b=s{,}'		=> \@barcodes,
	'hammingdist1|hd1=s'	=> \$hammingdist1,
	'hammingdist2|hd2=s'	=> \$hammingdist2,
	'help|h!'				=> \$help,
);

if($help){
print "
USAGE: $0
  -h: help
  -i: input data table with barcodes ID, required
  -o: prefix of output files, default 'tmpdata'
  -b: barcode files (order must be I5 I7 T5 T7 for each sample in samples), default 'I5 I7 T5 T7'
  -bf: path to barcode files, default '~/04.snATAC.pig/barcode.lib'
  -hd1: mismatch allowed for each sequence from reference dict, default '2'
  -hd2: distance allowed for the minimal and second sequence, default '2'
 
Note:
  Two output files: tmpdata.demx.barcode.barcodes and tmpdata.demx.barcode.samples are input files for 'get_fixed.fq.indel.pl'.
  
  Returned barcode format: T7I7T5I5
  **split samples
  **example of input table: 
'
LW-2w-Liver-2	T5-1:T5-8	T7-1:T7-12	N301-320	N401-496
LW-2w-Heart-2	T5-9:T5-16	T7-13:T7-24	N309-336+N349-352	N401-496
'
  
  ********************output for barcodes:
  1. mismatchs for each barcode T5/T7/I5/I7 (2, 2)
  2. insertions ......

";

}elsif($input){
	$output||="tmpdata";
	$barcode_folder||="~/04.snATAC.pig/barcode.lib";
	@barcodes = ("I5","I7","T5","T7") if !@barcodes;
	$hammingdist1||=2;
	$hammingdist2||=2;
	
	# barcode hash
	my %i5 = barcode_hash("$barcode_folder/$barcodes[0]", "I");
	my %i7 = barcode_hash("$barcode_folder/$barcodes[1]", "I");
	my %t5 = barcode_hash("$barcode_folder/$barcodes[2]", "T");
	my %t7 = barcode_hash("$barcode_folder/$barcodes[3]", "T");
	
	# process data table
	# generate sample/barcode dict
	open IN,$input or die;
	open OUT,">$output.demx.barcode.samples" or die;
	my %i5_2;
	my %i7_2;
	my %t5_2;
	my %t7_2;
	while(<IN>){
		chomp;
		my @arr = split;
		foreach my $T7 (T57_array($arr[2],"T7")){
			$t7_2{$t7{$T7}} = $T7;
			foreach my $I7 (I57_array($arr[4])){
				$i7_2{$i7{$I7}} = $I7;
				foreach my $T5 (T57_array($arr[1],"T5")){
					$t5_2{$t5{$T5}} = $T5;
					foreach my $I5 (I57_array($arr[3])){
						$i5_2{$i5{$I5}} = $I5;
						print OUT "$arr[0]\t$T7$I7$T5$I5\n";
					}
				}
			}
		}
	}
	close IN;
	close OUT;
	
	# generate good barcode dict
	system("rm -rf $output.demx.barcode.barcodes");
	open my $outbarcodes,">>","$output.demx.barcode.barcodes" or die;
	# T7
	get_miss_bar("t7", $hammingdist1, $hammingdist2, $outbarcodes, \%t7_2);
	get_insert_T("t7", 1, $hammingdist2, $outbarcodes, \%t7_2);
	# I7
	get_miss_bar("i7", $hammingdist1, $hammingdist2, $outbarcodes, \%i7_2);
	get_insert_I("i7", "C", 1, $hammingdist2, $outbarcodes, \%i7_2);
	get_delete_I("i7", 1, $hammingdist2, $outbarcodes, \%i7_2);
	# T5
	get_miss_bar("t5", $hammingdist1, $hammingdist2, $outbarcodes, \%t5_2);
	get_insert_T("t5", 1, $hammingdist2, $outbarcodes, \%t5_2);
	# I5
	get_miss_bar("i5", $hammingdist1, $hammingdist2, $outbarcodes, \%i5_2);
	get_insert_I("i5", "A", 1, $hammingdist2, $outbarcodes, \%i5_2);
	get_delete_I("i5", 1, $hammingdist2, $outbarcodes, \%i5_2);
	close $outbarcodes;
	
}else{
	exec("$0 -h")
}

# construct barcode hash in reference lib
sub barcode_hash{
	my ($bcfile, $TI) = @_;
	my %hash;
	open IN,"$bcfile" or die;
	while(<IN>){
		chomp;
		my @arr = split;
		#$arr[1] =~ s/^.// if $TI eq "T";
		$hash{$arr[0]} = $arr[1];
	}
	close IN;
	return %hash;
}

# genrate barcode array for T5/T7 of sample
sub T57_array{
	my ($t57_in, $T57) = @_;
	my @t57_arr;
	foreach my $group1 (split('\+',$t57_in)){
		my @group2 = split(":",$group1);
		if($#group2==0){
			push @t57_arr,$group2[0];
		}elsif($#group2==1){
			my $start = (split("-",$group2[0]))[1];
			my $end = (split("-",$group2[1]))[1];
			for($start..$end){
				push @t57_arr,"$T57-$_";
			}
		}else{
			die("Wrong barcode format in table: $group1\n");
		}
	}
	return @t57_arr;
}

# genrate barcode array for I5/I7 of sample
sub I57_array{
	my $i57_in = @_[0];
	my @i57_arr;
	foreach my $group1 (split('\+',$i57_in)){
		my @group2 = split("-",$group1);
		if($#group2==0){
			push @i57_arr,$group2[0];
		}elsif($#group2==1){
			$group2[0] =~ s/^N//;
			for($group2[0]..$group2[1]){
				push @i57_arr,"N$_";
			}
		}else{
			die("Wrong barcode format in table: $group1\n");
		}
	}
	return @i57_arr;
}

## define a subrotine hd for hamming distance
sub hd {
	return ($_[0] ^ $_[1]) =~ tr/\001-\255//;
}

# get mismatched barcodes
sub get_miss_bar{
	my ($name, $minimal, $second, $fh_out, $hash_barcode) = @_;
	foreach my $seqarr (cartesian {return @_} [qw(A T C G)] x 8){
		my $seq = join '',@$seqarr;  #convert array index to characters
		print $fh_out "$name\t$hash_barcode->{$seq}\t$seq\t$seq\t0\tM\n" and next if exists $hash_barcode->{$seq};
		my $mindis = 100;
		my $seconddis = 100;
		my $matchBD;
		foreach my $key (keys %{$hash_barcode}){
			my $hamdis = hd($seq,$key);
			if($hamdis < $mindis){
				$seconddis = $mindis;
				$mindis = $hamdis;
				$matchBD = $key;
			}elsif($hamdis >= $mindis and $hamdis < $seconddis){
				$seconddis = $hamdis;
			}
		}
		print $fh_out "$name\t$hash_barcode->{$matchBD}\t$seq\t$matchBD\t$mindis\tM\n" if $mindis <= $minimal and $seconddis-$mindis >= $second;
	}
}

# get insertion barcodes for T5/T7
sub get_insert_T{
	my ($name, $minimal, $second, $fh_out, $hash_barcode) = @_;
	## get new hash
	my %newbarcodes;
	my %newbarcodes2;
	foreach (keys %{$hash_barcode}){
		my $key = substr($_, 0, 7);
		$newbarcodes{$key} = $hash_barcode->{$_};
		$newbarcodes2{$key} = $_;
	}
	## generate barcodes
	my @bases = ("A", "G", "C", "T");
	foreach my $seqarr (cartesian {return @_} [qw(A T C G)] x 7){
		my $seq = join '',@$seqarr;
		if(exists $newbarcodes{$seq}){
			foreach my $base (@bases){
				print $fh_out "$name\t$newbarcodes{$seq}\t$base$seq\t$newbarcodes2{$seq}\t0\tI\n";
			}
			next;
		}
		my $mindis = 100;
		my $seconddis = 100;
		my $matchBD;
		foreach my $key (keys %newbarcodes){
			my $hamdis = hd($seq,$key);
			if($hamdis < $mindis){
				$seconddis = $mindis;
				$mindis = $hamdis;
				$matchBD = $key;
			}elsif($hamdis >= $mindis and $hamdis < $seconddis){
				$seconddis = $hamdis;
			}
		}
		if($mindis <= $minimal and $seconddis-$mindis >= $second){
			foreach my $base (@bases){
				print $fh_out "$name\t$newbarcodes{$matchBD}\t$base$seq\t$newbarcodes2{$matchBD}\t$mindis\tI\n";
			}
		}
	}
}

# get insertion barcodes for I5/I7
sub get_insert_I{
	my ($name, $fixindex_last, $minimal, $second, $fh_out, $hash_barcode) = @_;
	## get new hash
	my %newbarcodes;
	my %newbarcodes2;
	foreach (keys %{$hash_barcode}){
		my $key = $fixindex_last . substr($_, 0, 7);
		$newbarcodes{$key} = $hash_barcode->{$_};
		$newbarcodes2{$key} = $_;
	}
	## generate barcodes
	foreach my $seqarr (cartesian {return @_} [qw(A T C G)] x 8){
		my $seq = join '',@$seqarr;
		print $fh_out "$name\t$newbarcodes{$seq}\t$seq\t$newbarcodes2{$seq}\t0\tI\n" and next if exists $newbarcodes{$seq};
		my $mindis = 100;
		my $seconddis = 100;
		my $matchBD;
		foreach my $key (keys %newbarcodes){
			my $hamdis = hd($seq,$key);
			if($hamdis < $mindis){
				$seconddis = $mindis;
				$mindis = $hamdis;
				$matchBD = $key;
			}elsif($hamdis >= $mindis and $hamdis < $seconddis){
				$seconddis = $hamdis;
			}
		}
		print $fh_out "$name\t$newbarcodes{$matchBD}\t$seq\t$newbarcodes2{$matchBD}\t$mindis\tI\n" if $mindis <= $minimal and $seconddis-$mindis >= $second;
	}
}

# get deletion barcodes for I5/I7
sub get_delete_I{
	my ($name, $minimal, $second, $fh_out, $hash_barcode) = @_;
	## get new hash
	my %newbarcodes;
	my %newbarcodes2;
	foreach (keys %{$hash_barcode}){
		my $key = substr($_, 1, 7);
		$newbarcodes{$key} = $hash_barcode->{$_};
		$newbarcodes2{$key} = $_;
	}
	## generate barcodes
	my @bases = ("A", "G", "C", "T");
	foreach my $seqarr (cartesian {return @_} [qw(A T C G)] x 7){
		my $seq = join '',@$seqarr;
		if(exists $newbarcodes{$seq}){
			foreach my $base (@bases){
				print $fh_out "$name\t$newbarcodes{$seq}\t$seq$base\t$newbarcodes2{$seq}\t0\tD\n";
			}
			next;
		}
		my $mindis = 100;
		my $seconddis = 100;
		my $matchBD;
		foreach my $key (keys %newbarcodes){
			my $hamdis = hd($seq,$key);
			if($hamdis < $mindis){
				$seconddis = $mindis;
				$mindis = $hamdis;
				$matchBD = $key;
			}elsif($hamdis >= $mindis and $hamdis < $seconddis){
				$seconddis = $hamdis;
			}
		}
		if($mindis <= $minimal and $seconddis-$mindis >= $second){
			foreach my $base (@bases){
				print $fh_out "$name\t$newbarcodes{$matchBD}\t$seq$base\t$newbarcodes2{$matchBD}\t$mindis\tD\n";
			}
		}
	}
}

