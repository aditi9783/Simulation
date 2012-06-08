
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
my $snpset;
my $errthres = 0.2;			# min avg number of links for a position to be considered as a true hap position
my @newindex;
my @haps;
my @newpos;
my $numsnps;
my ($maxid, $newinfo, $conflict, $overlap);
my (@idx, @sorted);
my ($expidx, $hid, $numgaps);

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
#open (SET, ">$snpfile.snpset_noerrsnp" ) or die "can't open file to write: $!\n";
#open (SS, ">$snpfile.sorted_snpset_noerrsnp_reindexed" ) or die "can't open file to write: $!\n";
open (HAP, ">$snpfile.hapsfromsnpsets" ) or die "can't open file to write: $!\n";

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
	my ($pflag, $rid, $idx, @read);

	while ( <FH> ) {
		if ($_ =~ /\>p(\d+)\s\:r\s(\d+)\s\:idx\s(\d+)/) {
			$pflag = $1;	# pflag = 0 for first read and = 1 for paired read
			$rid = $2;
			$idx = $3;
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
			# inset a marker to signal end of a read (useful to distinguish paired reads)
			push @{ $snpreads->{$rid} }, "E" if ($pflag == 0);
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
			$index1 = $snpreads->{$i}->[$j];
			next if $index1 eq "E";				# marker element: signals read end
			for $k ( $j+1  .. $#{ $snpreads->{$i} } ) {
				$index2 = $snpreads->{$i}->[$k];
				next if $index2 eq "E";
				$snp_pairs[$index1][$index2]++;
				$snp_pairs[$index2][$index1]++;		# to make matrix symmetric
			}
		}
	}
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
		if ($snp_pos[$i+1] == $snp_pos[$i]) {	# if a position is mapped with different bases (i.e. position occurs more than once), then its likely not an error position.
			$i += 2;
		} else {
			push @errpos, $i if $avg <= $errthres;
		}
		# re-index the nonerror positions
		push @newindex, $i if $avg > $errthres;
#		print SUM "$i: $snp_pos[$i] $snp_base[$i] $sum\n" if $avg > $errthres;		# print number of links for non error positions
	}
	# free @snp_pairs
	undef(@snp_pairs);
	for $i ( keys %{ $snpreads } ) {
		@temp = ();
		for $j ( 0 .. $#{ $snpreads->{$i} } ) {
			# reindex non-error positions so that they are consecutively numbered
			push @temp, "E" if ( $snpreads->{$i}->[$j] eq "E" );
			for $k ( 0 .. $#newindex ) {
				if ( $snpreads->{$i}->[$j] == $newindex[$k] ) {
					push @temp, $k;
					last;
				}
			}
		}
		if ($#temp >= 2) {	# there are at least two non-error snps (and a marker "E") that are present in this read
#			$flag = 0;
			if ( $temp[0] eq "E" ) {	# snps present in paired read only (not in first read)
				# test if the current snpset is already present
#				for $k ( 0 .. $#{ $snpset->[ $temp[1] ] } ) {
#					if ( @temp[1 .. $#temp] ~~ @{ $snpset->[ $temp[1] ]->[$k] } ) {
#						$flag = 1;
#						last;
#					}
#				}
				push @{ $snpset->[ $temp[1] ] }, [ @temp[1 .. $#temp] ];
			} else {
				push @{ $snpset->[ $temp[0] ] }, [ @temp ];
			}
		}
	}
	undef ($snpreads);
	$hid = 0;	# hap id
	my (@matches, $maxhid); 
	for $i ( 0 .. $#{ $snpset } ) {
		next if $#{ $snpset->[$i] } < 0;	# no snpsets for this position, probably an error site
		for $j ( 0 .. $#{ $snpset->[$i] } ) {
			# store indexes for snps that are missing in snpsets
			@temp = ();
			$expidx = -1;
#			print HAP "@{ $snpset->[$i]->[$j] }\n" if $i == 0;
			for $k ( 0 .. $#{ $snpset->[$i]->[$j] } ) {
				if ( $snpset->[$i]->[$j]->[$k] eq "E" ) {
					$expidx = -1;
				}
				elsif ($expidx == -1 || $expidx == $snpset->[$i]->[$j]->[$k]) {
					$numsnps++;
					$expidx = $snpset->[$i]->[$j]->[$k]+1 ;		# consecutive indexes are expected, unless marker "E" is seen	
				} else {
					$numgaps = $snpset->[$i]->[$j]->[$k] - $expidx;	# number of consecutive snps that are missing
					while ( $numgaps > 0 ) {
						push @temp, $snpset->[$i]->[$j]->[$k] - $numgaps;	# snp is absent from this set
						$numgaps--;
					}
					$expidx = $snpset->[$i]->[$j]->[$k]+1 ;		# consecutive indexes are expected
				}
			}
			next if $#temp < 0;	# no gaps in this snpset
			print "$i: @temp\n";
			next;
			# when comparing @temp to exisitng haps: three possibilities: conflict, expands current hap, is contained in current hap
			# if conflict: add a new row, else extend an exisiting hap
			push @haps, [ @temp ] if $#haps < 0;	# insert first snp set
			
			# check if current hap can be merged with existing ones, or if there is a conflict
			updateHap( \@temp );
		}
	}
	# print index and corresponding snp details
	for $i ( 0 .. $#newindex ) {
		# old index is $newindex[$i]
		print HAP "$i: $snp_pos[ $newindex[$i] ] $snp_base[ $newindex[$i] ]\n";
	}
	print HAP "Number of haps: $#haps+1\n";
	for $i ( 0 .. $#haps ) {
		print HAP "$i:\t @{ $haps[$i] }\n\n";
	}
	$snp_pos[ $#snp_pos+1 ] = undef;	# remove dummy var from array	
}
# ENd of findErr #

##################
# updateHap
#
# takes in a list of snp indexes that are absent in a read, and see
# if they can be merged with any exisinting hap deletion profiles
##################
sub updateHap {
	my $del = shift;

	for my $a ( 0 .. $#haps ) {
		for my $b ( 0 .. $#{ $haps->[$a] } ) {
			for my $c ( 0 .. $#del ) {
				
			}	
		}
	}
}
# End of updateHap #
