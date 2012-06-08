#!/usr/bin/perl;

use strict;
use Getopt::Std;

our ($opt_r, $opt_e, $opt_t, $opt_h);
getopt( 'reth' );

my $t = 3;
my $haps;
my @unique_hpos;	# unique hap positions in all the haplotypes
my @h;

$t = $opt_t if $opt_t;

if ($opt_r and $opt_e) {
	open(FH, $opt_r ) or die "can't open file $opt_r: $!\n";
	open( ER, $opt_e ) or die "can't open file $opt_e: $!\n";
	open( OUT, ">$opt_e.$t.poissoncorr" ) or die "can't open file to write: $!\n";
	open( SEL, ">$opt_e.$t.selectPos" ) or die "can't open file to write: $!\n";
	open( SNP, ">$opt_e.$t.snp" ) or die "can't open file to write: $!\n";
} else {
	print "USAGE\n-r\tReference file in fasta format (required).\n";
	print "-e\tError file containing error bases at each position of mapped reads (required).\n";
	print "-t\tThreshold for pruning error bases at each pos (default 3). Threshold is how many SDs from mean should the count be to be selected.\n";
	print "-h\tHaplotype file to calculate precision and recall of positions recovered.\n";
	die "\n"; 
}

# to calculate precision and recall of positions recovered as opposed to true haplotype positions
if ($opt_h) {
	open( HAP, $opt_h ) or die "can't open file $opt_h: $!\n";
	while ( <HAP> ) {
		if ($_ =~ /.+\t\d+\t[A-Z]/) {
			@h = split /\t/;
			for (my $i=1; $i < $#h; $i +=2 ) {
				push @{ $haps }, [ $h[0], $h[$i] ];	# h[0] contains hap freq
				push @unique_hpos, $h[$i];
			}
		}
	}
}

#remove redundant hap positions
my @sorted_upos = sort { $a <=> $b } @unique_hpos;
@unique_hpos = ();
push @unique_hpos, $sorted_upos[0];
print "Unique hap positions:\n@unique_hpos";
for my $i ( 1 .. $#sorted_upos ) {
	next if ( $sorted_upos[$i] == $sorted_upos[$i-1] );
	push @unique_hpos, $sorted_upos[$i];
	print "$sorted_upos[$i]\n";
}


my @selectPos; 		# selected positions after poisson correction
my $seq;
while (<FH> ) {
	if ($_ =~ /^[a-zA-Z]+/) {
		chomp;
		$seq = $seq . $_;
	}
}

my @genome = split "", $seq;

my @base = ("A", "T", "C", "G");
my (@count, @poisson, @orig_count);
my ($nA, $nT, $nC, $nG) = (0,0,0,0);
my ($sum, $corrsum, $thres, $exp);

#print OUT "pos: truebase - orig sum - orig count - exp count - corr sum - corr count - poisson prob\n";
while (<ER>) {
	if ($_ =~ /(\d+)\:\t0\s(.+)\n/ ) {
	#	print OUT "\n";
		my @a = split /\s/, $2;
		my $b = join "", @a;
		$nA = ($b =~ tr/A/A/ );
		$nT = ($b =~ tr/T/T/ );
		$nC = ($b =~ tr/C/C/ );
		$nG = ($b =~ tr/G/G/ );
		@count = ($nA, $nT, $nC, $nG);
		@orig_count = @count;
		$sum = 0;

		for my $j ( 0 .. $#base ) {
			$count[ $j ] = 0 if $base[$j] eq $genome[$1 - 1];	# ignore error bases that match the true base
			$sum += $count[$j];
		}
		print OUT "$1: $genome[$1 - 1]\t$sum\t@count";
		next if $sum == 0;
		
	#	poissonCal ();
		$corrsum = 0;

		$exp = $sum/3;	# expecting uniform distribution for the three errored bases
		# $exp is also the mean of the normal approx of the poisson for counts > 100. SD = sqrt( $exp ).
		# identify those bases whose occurence is 2 SD greater than mean as putative hap positions

		$thres = $exp + ($t * sqrt( $exp ));

		for my $j ( 0 .. $#base ) {
			if ( $count[$j] < $thres ) {
				$count[$j] = 0;		# not a hap base
			} else {
				print SNP "$1\t$base[$j]\t$count[$j]\n";
			}
			$corrsum += $count[$j];
		}
		print OUT "-\t$corrsum - @count\n";
		print SEL "$1: @orig_count - $corrsum - @count\n" if $corrsum > 0;
		push @selectPos, $1 if $corrsum > 0;
	}
}

close ER;

# calculate precision and recall
my ($tp, $fp) = (0,0);
my $flag = 0;
my $nf_infreq = 0;	# number of hap positions not found that have freq < 0.001
if ($opt_h) {
	# calculate unique hap pos that are identified
	for my $upos ( @unique_hpos ) {
		for my $pos ( @selectPos ) {
			$tp++ if $pos == $upos;
		}
	}
	
	# find hap pos and their frequncies that are not identified.
	for my $i ( 0 .. $#{ $haps } ) {
		$flag = 0;
		for my $pos ( @selectPos ) {
			if ($pos == $haps->[$i]->[1]) {
				$flag = 1;
				last;
			}
		}
		print "NOT FOUND: @{ $haps->[$i] }\n" if $flag == 0;
		if ($flag == 0 and $haps->[$i]->[0] >= 0.001) {
			$nf_infreq++;
		}
	}
}

my $precision = $tp/($#selectPos+1);
#my $recall = $tp/($#{ $haps } +1);
my $recall = $tp/($#unique_hpos+1);

print "true hap positions identified: $tp\nNum positions selected: $#selectPos+1\tTotal num of (unique) true hap pos: $#unique_hpos+1\n";
print "Precision: $precision\tRecall: $recall\n";
print "Number of hap positions that have freq >= 0.001 and were not found: $nf_infreq\n";
################

sub factorial {
  my $n = shift;
  $n == 0 ? 1 : $n*factorial($n-1);
}

################

sub poissonCal {
		@poisson = (0,0,0,0);
		$thres = $t/$sum;
		for my $j ( 0 .. $#base ) {
			next if $base[$j] eq $genome[$1 - 1];
			$poisson[$j] = ($exp ** $count[$j]) * exp( -1 * $exp ) / factorial( $count[$j] );
			if ( $poisson[$j] > $thres ) {   
				$count[$j] = 0;
			}
			$count[$j] = 0 if ( $exp > $count[$j] );	# remove base occurences that are rarer than the expected
			$corrsum += $count[$j];
		}
		print OUT "\t$corrsum\t@count\t\t@poisson";
}	
