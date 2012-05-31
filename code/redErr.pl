#!/usr/bin/perl;

use strict;
use Getopt::Std;

our ($opt_r, $opt_e, $opt_t, $opt_h, $opt_c);
getopt( 'rethc' );

my @base = ("A", "T", "C", "G");
my ($nA, $nT, $nC, $nG) = (0,0,0,0);
my (%snp, %pred_snp);
my (@genome, @coverage);

if ($opt_r and $opt_e) {
	open( ER, $opt_e ) or die "can't open file $opt_e: $!\n";
	open( OUT, ">$opt_e.poissoncorr" ) or die "can't open file to write: $!\n";
	open( SNP, ">$opt_e.snp" ) or die "can't open file to write: $!\n";
} else {
	print "USAGE\n-r\tReference file in fasta format (required).\n";
	print "-e\tError file containing error bases at each position of mapped reads (required).\n";
	print "-h\tHaplotype file to calculate precision and recall of positions recovered.\n";
	print "-c\tCoverage file (required).\n";
	die "\n"; 
}

# to calculate precision and recall of positions recovered as opposed to true haplotype positions
getSNPs() if ($opt_h);
getSeq();	# get reference seq from file and store in @genome
getCoverage(); 	# get coverage for each position from file and store it in @coverage

while (<ER>) {
	if ($_ =~ /(\d+)\:\t0\s(.+)\n/ ) {
		my $pos = $1;
		my @a = split /\s/, $2;
		my $b = join "", @a;
		$nA = ($b =~ tr/A/A/ );
		$nT = ($b =~ tr/T/T/ );
		$nC = ($b =~ tr/C/C/ );
		$nG = ($b =~ tr/G/G/ );
		my @count = ($nA, $nT, $nC, $nG);
		my ($sum, $corrsum) = (0,0);

		for my $j ( 0 .. $#base ) {
			$count[ $j ] = 0 if $base[$j] eq $genome[$pos - 1];	# ignore error bases that match the true base
			$sum += $count[$j];
		}
		print OUT "$pos: $genome[$pos - 1]\t$sum\t@count";
		next if $sum == 0;

		my $exp = $sum/3;	# expecting uniform distribution for the three errored bases
		# $exp is also the mean of the normal approx of the poisson for counts > 100. SD = sqrt( $exp ).
	
		my $t = 0.012 * sqrt($coverage[$pos]);		# sd that will identify a base that is a snp at freq 0.001
		my $thres = $exp + ($t * sqrt( $exp ));	# true snps should have counts > threshold
		
		my $minbasecount = 0.00385 * $coverage[$pos]; # a pos with snp at 0.001 theoretically has base count of 0.0043 (not considering errors in haps)
							      # but adjusted to 0.00385 heuristically
		my $thresflag = 0;

		for my $j ( 0 .. $#base ) {
			if ( $count[$j] >= $thres ) {
				print SNP "$pos\t$base[$j]\t$count[$j]\t$coverage[$pos]\n";
				$pred_snp{$pos.$base[$j]}++;
				$thresflag = 1;				# snp is identified by threshold method
				$corrsum += $count[$j];
			}
		}
		if ($thresflag == 0) {
			for my $j ( 0 .. $#base ) {
				next if $count[$j] < $minbasecount;
				print SNP "$pos\t$base[$j]\t$count[$j]\t$coverage[$pos]\n";
				$pred_snp{$pos.$base[$j]}++;
				$corrsum += $count[$j];
			}		
		}
		print OUT "-\t$corrsum\n";
	}
}
close ER;

# calculate precision and recall for acutal snps, i.e. snps at same positions are considered separately
my $correctpred = 0;	# number of true snps identified as true snps
my $falsesnp = 0;	# number of true snps missed
my $totalsnp = 0;	# total number of true snp pos
my $totalpred = 0;	# total number of sites predicted

for my $s ( keys %snp ) {
	$totalsnp++;
	print "\ntrue snp: $s: $snp{$s}\t";
	if (defined($pred_snp{$s})) {
		$correctpred++;
	} else {
		print "-- missed\n";
	}	
}

for my $ps ( keys %pred_snp ) {
	$totalpred += $pred_snp{$ps};
}

$falsesnp = $totalsnp - $correctpred;
print "\n\nFrom actual snps:\nTotal snps: $totalsnp\nTotal predicted snps: $totalpred\nCorrectly predicted snps: $correctpred\nMissed snps: $falsesnp\n";

my $precision = $correctpred/$totalpred;
my $recall = $correctpred/$totalsnp;
print "\nPrecision (tp/total pred) = $precision\n";
print "Recall (tp/total snp) = $recall\n";

#############
# getSNPs
#
# get SNPs from the hap file
#############
sub getSNPs {
	open( HAP, $opt_h ) or die "can't open file $opt_h: $!\n";
	while ( <HAP> ) {
		if ($_ =~ /.+\t\d+\t[A-Z]/) {
			chomp;
			my @h = split /\t/;
			for (my $i=1; $i < $#h; $i +=2 ) {
				$snp{$h[$i].$h[$i+1]}++;	# snp pos and base
			}
		}
	}
	close HAP;
}
# End of getSNPs #

##############
# readSeq
#
# read genome seq from a fasta file
##############
sub getSeq {
	open(FH, $opt_r ) or die "can't open file $opt_r: $!\n";
	my $seq;
	while (<FH> ) {
		if ($_ =~ /^[a-zA-Z]+/) {
			chomp;
			$seq = $seq . $_;
		}
	}
	close FH;
	@genome = split "", $seq;	# @genome is a global variable that now stores the seq
}
# End of getSeq #

################
# getCoverage
###############
sub getCoverage {
	open( COV, $opt_c ) or die "can't open coverage file $opt_c: $!\n";
	$coverage[0] = "dummy";		# for coverage at a seq pos, index starts from 1
	while ( <COV> ) {
		if ($_ =~ /(\d+)\:\t(\d+)/ ) {
			$coverage[ $1 ] = $2;
		}
	}
	close COV;
}
# End of getCoverage #
