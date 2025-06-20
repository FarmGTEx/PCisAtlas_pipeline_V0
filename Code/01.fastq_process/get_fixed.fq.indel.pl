#!/usr/bin/env perl
#pre-processing fastq files of snATAC data
use Getopt::Long;
use Array::Utils qw(:all);
use List::MoreUtils ':all';

GetOptions(
	'inputs|i=s{,}'			=> \@inputs,
	'barcodes1|b1=s'			=> \$barcodes1,
	'barcodes2|b2=s'			=> \$barcodes2,
	'output|o=s'			=> \$output,
	'thread|t=i'			=> \$thread,
	'help|h!'			=> \$help,
);

if($help){
print "\n";
print "usage: $0\n";
print " -h: help\n";
print " -i: input fastq files [I1, I2, R1, R2]\n";
print " -b1: barcode file 1\n";
print " -b2: barcode file 2\n";
print " -t: threads used in compression step, default '1'\n";
print " -o: prefix of output fastqs, default 'tmpdata'\n";

print "
Note:
  Returned barcode format: T7I7T5I5
  **fixed index are not considered in this script, only the four pieces of barcodes are included

";

}elsif(@inputs and $barcodes1 and $barcodes2){
	$output||="tmpdata";
	$thread||=1;
	
	##barcode hash
	open IN,$barcodes1 or die;
	my (%hash_i5_m,%hash_i7_m,%hash_t5_m,%hash_t7_m,%hash_i5_i,%hash_i7_i,%hash_t5_i,%hash_t7_i,%hash_i5_d,%hash_i7_d) = ('','','','','','','','','','');
	while(<IN>){
		chomp;
		my @arr = split;
		if($arr[0] eq "i5" and $arr[5] eq "M"){
			$hash_i5_m{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "i5" and $arr[5] eq "I"){
			$hash_i5_i{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "i5" and $arr[5] eq "D"){
			$hash_i5_d{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "i7" and $arr[5] eq "M"){
			$hash_i7_m{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "i7" and $arr[5] eq "I"){
			$hash_i7_i{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "i7" and $arr[5] eq "D"){
			$hash_i7_d{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "t5" and $arr[5] eq "M"){
			$hash_t5_m{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "t5" and $arr[5] eq "I"){
			$hash_t5_i{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "t7" and $arr[5] eq "M"){
			$hash_t7_m{$arr[2]} = [$arr[1],$arr[3]];
		}elsif($arr[0] eq "t7" and $arr[5] eq "I"){
			$hash_t7_i{$arr[2]} = [$arr[1],$arr[3]];
		}else{
			die("Wrong format of barcode dict file!");
		}
	}
	close IN;
	# sample dict
	my %samples;
	open IN,$barcodes2 or die;
	while(<IN>){
		chomp;
		my @arr = split;
		$samples{$arr[1]} = $arr[0];
	}
	close IN;
	
	##extracting good fastq files
	if($inputs[0] =~ /gz$/){
		open INDEX1,"unpigz -dkc -p $thread $inputs[0]|" or die;
		open INDEX2,"unpigz -dkc -p $thread $inputs[1]|" or die;
		open FQ1,"unpigz -dkc -p $thread $inputs[2]|" or die;
		open FQ2,"unpigz -dkc -p $thread $inputs[3]|" or die;
	}else{
		open INDEX1,"$inputs[0]" or die;
		open INDEX2,"$inputs[1]" or die;
		open FQ1,"$inputs[2]" or die;
		open FQ2,"$inputs[3]" or die;
	}
	my %OUTFH;
	foreach my $sample (uniq values %samples){
		open $OUTFH{"$sample.r1"},">${output}.demulti.$sample\_R1.fastq" or die;
		open $OUTFH{"$sample.r2"},">${output}.demulti.$sample\_R2.fastq" or die;
	}
	my $allreads;
	my $goodall;
	my @stat_gorup_reads;  #T5, I5, T7, I7
	my %sample_count;
	while(!eof(INDEX1) and !eof(INDEX2) and !eof(FQ1) and !eof(FQ2)){
		$allreads++;
		#INDEX1
		my $I1line1=<INDEX1>;
		my $I1line2=<INDEX1>;
		my $I1line3=<INDEX1>;
		my $I1line4=<INDEX1>;
		#INDEX2
		my $I2line1=<INDEX2>;
		my $I2line2=<INDEX2>;
		my $I2line3=<INDEX2>;
		my $I2line4=<INDEX2>;
		#FQ1
		my $F1line1=<FQ1>;
		my $F1line2=<FQ1>;
		my $F1line3=<FQ1>;
		my $F1line4=<FQ1>;
		#FQ2
		my $F2line1=<FQ2>;
		my $F2line2=<FQ2>;
		my $F2line3=<FQ2>;
		my $F2line4=<FQ2>;
		
		# assign to correct cell barcode (barcode correction)
		my $T7 = substr $I1line2,0,8;
		my $I7 = substr $I1line2,35,8;
		my $T5 = substr $I2line2,0,8;
		my $I5 = substr $I2line2,29,8;
		my ($T7_correct, $I7_correct, $T5_correct, $I5_correct) = ('','','','');
		my ($T7_id, $I7_id, $T5_id, $I5_id) = ('','','','');
		if(exists $hash_i5_m{$I5} and exists $hash_i7_m{$I7} and exists $hash_t5_m{$T5} and exists $hash_t7_m{$T7}){
			$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_m{$I7}[1]; $T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_m{$I5}[1];
			$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_m{$I7}[0]; $T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_m{$I5}[0];
			$stat_gorup_reads[0] += 1;
		}elsif(exists $hash_i5_m{$I5} and exists $hash_t5_m{$T5}){
			$T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_m{$I5}[1];
			$T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_m{$I5}[0];
			if($hash_t7_m{$T7} and exists $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[1]++;
			}elsif(exists $hash_t7_m{$T7} and $hash_i7_d{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_d{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_d{$I7}[0];
				$stat_gorup_reads[2]++;
			}elsif(exists $hash_t7_i{$T7} and $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_i{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_i{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[3]++;
			}else{
				$bad_reads++;
				next;
			}
		}elsif(exists $hash_i7_m{$I7} and exists $hash_t7_m{$T7}){
			$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_m{$I7}[1];
			$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_m{$I7}[0];
			if($hash_t5_m{$T5} and exists $hash_i5_i{$I5}){
				$T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_i{$I5}[1];
				$T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_i{$I5}[0];
				$stat_gorup_reads[4]++;
			}elsif(exists $hash_t5_m{$T5} and $hash_i5_d{$I5}){
				$T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_d{$I5}[1];
				$T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_d{$I5}[0];
				$stat_gorup_reads[8]++;
			}elsif(exists $hash_t5_i{$T5} and $hash_i5_i{$I5}){
				$T5_correct = $hash_t5_i{$T5}[1]; $I5_correct = $hash_i5_i{$I5}[1];
				$T5_id = $hash_t5_i{$T5}[0]; $I5_id = $hash_i5_i{$I5}[0];
				$stat_gorup_reads[12]++;
			}else{
				$bad_reads++;
				next;
			}
		}elsif(exists $hash_t5_m{$T5} and exists $hash_i5_i{$I5}){
			$T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_i{$I5}[1];
			$T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_i{$I5}[0];
			if($hash_t7_m{$T7} and exists $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[5]++;
			}elsif(exists $hash_t7_m{$T7} and $hash_i7_d{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_d{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_d{$I7}[0];
				$stat_gorup_reads[6]++;
			}elsif(exists $hash_t7_i{$T7} and $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_i{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_i{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[7]++;
			}else{
				$bad_reads++;
				next;
			}
		}elsif(exists $hash_t5_m{$T5} and exists $hash_i5_d{$I5}){
			$T5_correct = $hash_t5_m{$T5}[1]; $I5_correct = $hash_i5_d{$I5}[1];
			$T5_id = $hash_t5_m{$T5}[0]; $I5_id = $hash_i5_d{$I5}[0];
			if($hash_t7_m{$T7} and exists $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[9]++;
			}elsif(exists $hash_t7_m{$T7} and $hash_i7_d{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_d{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_d{$I7}[0];
				$stat_gorup_reads[10]++;
			}elsif(exists $hash_t7_i{$T7} and $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_i{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_i{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[11]++;
			}else{
				$bad_reads++;
				next;
			}
		}elsif(exists $hash_t5_i{$T5} and exists $hash_i5_i{$I5}){
			$T5_correct = $hash_t5_i{$T5}[1]; $I5_correct = $hash_i5_i{$I5}[1];
			$T5_id = $hash_t5_i{$T5}[0]; $I5_id = $hash_i5_i{$I5}[0];
			if($hash_t7_m{$T7} and exists $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[13]++;
			}elsif(exists $hash_t7_m{$T7} and $hash_i7_d{$I7}){
				$T7_correct = $hash_t7_m{$T7}[1]; $I7_correct = $hash_i7_d{$I7}[1];
				$T7_id = $hash_t7_m{$T7}[0]; $I7_id = $hash_i7_d{$I7}[0];
				$stat_gorup_reads[14]++;
			}elsif(exists $hash_t7_i{$T7} and $hash_i7_i{$I7}){
				$T7_correct = $hash_t7_i{$T7}[1]; $I7_correct = $hash_i7_i{$I7}[1];
				$T7_id = $hash_t7_i{$T7}[0]; $I7_id = $hash_i7_i{$I7}[0];
				$stat_gorup_reads[15]++;
			}else{
				$bad_reads++;
				next;
			}
		}else{
			$bad_reads++;
			next;
		}
		
		# get corrected barcodes and id
		my $thisbarcode = $T7_correct . $I7_correct . $T5_correct . $I5_correct;
		my $thisbarcode_id = $T7_id . $I7_id . $T5_id . $I5_id;
		next if not exists $samples{$thisbarcode_id};
		
		# output
		$F1line1 =~ s/^@//;
		$F2line1 =~ s/^@//;
		print {$OUTFH{"$samples{$thisbarcode_id}.r1"}} "\@$thisbarcode:$F1line1$F1line2$F1line3$F1line4";
		print {$OUTFH{"$samples{$thisbarcode_id}.r2"}} "\@$thisbarcode:$F2line1$F2line2$F2line3$F2line4";
		$goodall++;
		$sample_count{$samples{$thisbarcode_id}}++;
	}
	
	## output stat
	my $ratio_b = $goodall/$allreads*100;
	open OUT,">$output.fqstat" or die;
	print OUT "All read pairs: $allreads
Read pairs with corrected cell barcodes: $goodall ($ratio_b%)\nRead pairs of each sample:\n";
	foreach my $key (keys %sample_count){
		print OUT "  $key\t$sample_count{$key}\n";
	}
	print OUT "\nReads in each condition:\n";
	print OUT join("\n", @stat_gorup_reads) . "\n";  # order
	close OUT;
	
	## compress fastq files
	foreach my $key (keys %sample_count){
		system("pigz -p $thread ${output}.demulti.$key\_R1.fastq");
		system("pigz -p $thread ${output}.demulti.$key\_R2.fastq");
	}
	
}else{
	exec("$0 -h")
}

