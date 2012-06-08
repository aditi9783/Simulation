
#!/usr/bin/perl -w

# Reads in .snp file that has snp pos, snp base and # times it appears in reads (tab-separated). Looks up these snps in actual reads and finds out number of times any two snps are linked by same/paired read.

use strict;
use Getopt::Std;

our ($opt_a, $opt_b, $opt_s);
getopt( 'abs' );

my ($rfile1, $rfile2, $snpfile);	# read files in fasta format (rfile 1 and 2 for paired end reads); and snp file.
my (@snp_pos, @snp_base, @snp_count, @snp_pairs, $snpreads, @red_snppairs, @errpos);	# snps and their counts. 
										#the index of @snp will serve at indices for 2d array @snp_pairs to get mapping counts for any two snps.
										# @red_snppairs contains map counts for only those positions with avglinks > 0.2 or non-unique (when a position is 											part of > 1 haps using different bases, then one of the bases can have low #links).
										# hashref $snpreads has read id as key, and array of snps found in the read id as it value
										# @errpos has indices for snps that are likely errors and don't belong to true haps
my ($i, $j, $k, $index1, $index2, @temp, $sum, $avg);				# temporary variables
my ($flag, $set);
my %snpset;
my $errthres = 0.2;			# min avg number of links for a position to be considered as a true hap position

if ($opt_a and $opt_s) {
	$rfile1 = $opt_a;
	$snpfile = $opt_s;
} else {
	print "USAGE\n-a\tNGS read file in fasta format. Fasta header contains read id, paired or not info, and read start position. Required.\n";
	print "-b\tPaired NGS read file (paired reads for the reads provided with -a tag). Optional.\n";
	print "-s\tSNP file, contains snp position, snp base and #times snp occurs in all reads. Required\n";
	die "\n";
}

$rfile2 = $opt_b if $opt_b;

## MAIN ##

# get snps and their counts
getSNP();

# get read ids that have the snps
readSNP( $rfile1 );
readSNP( $rfile2 ) if $opt_b;

# initialize @snp_pairs
for $i ( 0 .. $#snp_pos ) {
	for $j ( 0 .. $#snp_pos ) {
		$snp_pairs[$i][$j] = 0;
	}
}

#open file handlers to store # links between snps (link matrix) and avg number of links per snp
#open (LNK, ">$snpfile.links" ) or die "can't open file to write: $!\n";
#open (SUM, ">$snpfile.numlink" ) or die "can't open file to write: $!\n";
open (SET, ">$snpfile.snpset_noerrsnp" ) or die "can't open file to write: $!\n";

# link snps based on shared read ids
linkSNP();

# find error positions based on how many times a position is linked with other positions
findErr();

# End of MAIN ##

#########################
# getSNP
#
# reads in snp file and returns snps
# in @snp and their counts in @snp_count
#########################
sub getSNP {
	open( SNP, $snpfile ) or die "can't open file $snpfile: $!\n";
	while ( <SNP> ) {
		chomp;
		@temp = split /\t/;
		push @snp_pos, $temp[0];
		push @snp_base, $temp[1];
		push @snp_count, $temp[2];	
	}	
	close SNP;
}
# End of getSNP #

#########################
# readSNP
#
# read file containing NGS reads in fasta format
# fasta header has read id and read start index.
# Returns hash reference where keys are read ids and 
# values are array ref containing snps in that read
#########################
sub readSNP {
	my $rfile = shift;
	open( FH, $rfile ) or die "can't open file $rfile: $!\n";
	my ($rid, $idx, @read);

	while ( <FH> ) {
		if ($_ =~ /\>.+r\s(\d+)\s\:idx\s(\d+)/) {
			$rid = $1;
			$idx = $2;
#			last if $rid > 5;
		} else {
			chomp;
			@read = split //;
			next if $#read < 0;	# the read is empty

			#check which snps are present in this read. SNP positions are sorted.
			for $i ( 0 .. $#snp_pos ) {
				last if $snp_pos[$i] > $idx+$#read;	# start index + readlen is the max sequence pos that the read can have
				next if $snp_pos[$i] < $idx;		# snp pos is less than read start pos
	
				push @{ $snpreads->{$rid} }, $i if ($read[ $snp_pos[$i] - $idx ] eq $snp_base[$i]);		
			}
		}
	}
	close FH;
}
# End of readSNP #

########################
# linkSNP
#
# reads the $snpreads that has read ids and snps therein
# and generates an adjaency matrix where the row and cols are 
# snp pos-snp base, and values indicate number of times a given
# pair of snps is found in same/paired read.
########################
sub linkSNP {
#	my @samereadsnp_hist;
	for $i ( keys %{ $snpreads } ) {
#		$samereadsnp_hist[ $#{ $snpreads->{$i} }+1 ]++;		# number of snps that share this read.
		for $j ( 0 .. $#{ $snpreads->{$i} }-1 ) {
			for $k ( $j+1  .. $#{ $snpreads->{$i} } ) {
				$index1 = $snpreads->{$i}->[$j];
				$index2 = $snpreads->{$i}->[$k];
				$snp_pairs[$index1][$index2]++;
				$snp_pairs[$index2][$index1]++;		# to make matrix symmetric
			}
		}
	}

#	print "Number of snps that are present in same read: Histogram\n";
#	for $i ( 0 .. $#samereadsnp_hist ) {
#		print "$i: $samereadsnp_hist[$i]\n" if $samereadsnp_hist[$i];
#	}

#	for my $set ( keys %snpset ) {
#		print SET "$snpset{ $set }: $set\n";	# number of times a set appears and the snp-set itself
#	}
	#print the adjacency matrix
#	print LNK "row and column indices:\n";
#	for $i ( 0 .. $#snp_pos ) {
#		print LNK "$i: $snp_pos[$i] $snp_base[$i]\n";
#	}
#	print LNK "\nSNP link matrix:\n";
#	for $i ( 0 .. $#snp_pos ) {
#		print LNK "@{ $snp_pairs[$i] }\n"; 
#	}
}

# End of linkSNP #

#######################
# findErr
#
# calculate number of times a position is linked with other positions
# (sum of each row in $snp_pairs), and returns the average number of links
# made overall. Error position should have fewer average # of links than
# real hap positions
#######################
sub findErr {
	$snp_pos[ $#snp_pos+1 ] = "X";		# dummy variable
	for $i ( 0 .. $#snp_pairs ) {
		($avg, $sum, $flag) = (0, 0, 0);
		for $j ( 0 .. $#{ $snp_pairs[$i] } ) {
			$sum += $snp_pairs[$i][$j];
		}
		$avg = $sum/( $#{$snp_pairs[$i]}+1 ) if $#{$snp_pairs[$i]} >=0 ;
#		print AVG "$snp_pos[$i]: $avg\n";
		
		if ($snp_pos[$i+1] == $snp_pos[$i]) {	# if a position is mapped with different bases (i.e. position occurs more than once), then its likely not an error position.
			$i += 2;
		} else {
			push @errpos, $i if $avg <= $errthres;
		}
#		print SUM "$i: $snp_pos[$i] $snp_base[$i] $sum\n" if $avg > $errthres;		# print number of links for non error positions
	}
	# free @snp_pairs
	undef(@snp_pairs);
	for $i ( keys %{ $snpreads } ) {
		@temp = ();
		for $j ( 0 .. $#{ $snpreads->{$i} } ) {
			# select non-error positions to be in the set
			$flag = 1;
			for my $err ( @errpos ) {
				$flag = 2 if $err == $snpreads->{$i}->[$j];
			}
			push @temp, $snpreads->{$i}->[$j] if $flag == 1;
		}
		if ($#temp >= 1) {	# there are at least two non-error snps that are present in this read
			$set = join( "-", @temp );
			$snpset{ $set }++;
		}
	}
	for $set (keys %snpset) {
		print SET "$snpset{ $set }: $set\n";
	}
	$snp_pos[ $#snp_pos+1 ] = undef;	# remove dummy var from array	
}
# ENd of findErr #
