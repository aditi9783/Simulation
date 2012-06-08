
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
	undef( $snpreads );
	$hid = 0;	# hap id
	my (@matches, $maxhid);
	for $i ( 0 .. $#{ $snpset } ) {
		next if $#{ $snpset->[$i] } < 0;	# no snpsets for this position, probably an error site
		for $j ( 0 .. $#{ $snpset->[$i] } ) {
			@temp = ();	# snps in this hap set, with missing snps marked by 'x'
			# intialize @temp with zeros, length = number of error-free snps
			for $k ( 0 .. $#newindex ) {
				$temp[$k] = 0;
			}
			$expidx = -1;
			$numsnps = 0;
#			print HAP "@{ $snpset->[$i]->[$j] }\n" if $i == 0;
			for $k ( 0 .. $#{ $snpset->[$i]->[$j] } ) {
				if ( $snpset->[$i]->[$j]->[$k] eq "E" ) {
					$expidx = -1;
					next;
				}
				elsif ($expidx == -1 || $expidx == $snpset->[$i]->[$j]->[$k]) {
					$temp[ $snpset->[$i]->[$j]->[$k] ] = 1;	
					$numsnps++;
					$expidx = $snpset->[$i]->[$j]->[$k]+1 ;		# consecutive indexes are expected, unless marker "E" is seen	
				} else {	# insert 'x' for gaps in snp set
					$numgaps = $snpset->[$i]->[$j]->[$k] - $expidx;	# number of consecutive snps that are missing
					while ( $numgaps > 0 ) {
						$temp[ $snpset->[$i]->[$j]->[$k] - $numgaps ] = "x";	# snp is absent from this set
						$numgaps--;
						$numsnps++;
					}
					$temp[ $snpset->[$i]->[$j]->[$k] ] = 1;	
					$numsnps++;
					$expidx = $snpset->[$i]->[$j]->[$k]+1 ;		# consecutive indexes are expected
				}
			}
			# when comparing @temp to exisitng haps: three possibilities: conflict, expands current hap, is contained in current hap
			# if conflict: add a new row, else extend an exisiting hap
			if ( $#haps < 0 ) {
				push @haps, [ @temp ];	# insert first snp set
#				print HAP "first hap:\n@temp\n" if $i == 0;
				next;			
			} else {
				@matches = ();	# stores number of snps that an exisitng hap shares with the new snpset. conflicted sets have matches = -1
				@newpos = ();
				for $hid ( 0 .. $#haps ) {
					for $k ( 0 .. $#newindex ) {
						if ( ($haps[$hid][$k] > 0 and $temp[$k] eq "x")	|| ($haps[$hid][$k] eq "x" and $temp[$k] > 0) ) {	# conflict
							$matches[$hid] = -1;
							last;
						} elsif ( ($haps[$hid][$k] > 0 and $temp[$k] > 0) || ($haps[$hid][$k] eq "x" and $temp[$k] eq "x") ) {
							$matches[$hid]++;
						} elsif ( ($haps[$hid][$k] == 0 and $temp[$k] > 0) || ($haps[$hid][$k] == 0 and $temp[$k] eq "x") ) {
							$newpos[$hid]++;	# new information this snpset has
						}	
					}
				}
				# extend the hap that has maximum matches and max newpos with the current snpset. If no matches, add snpset as a new hap (row)
				$maxhid = max( \@matches, \@newpos, $numsnps ); 
#				print HAP "new information wrt haps: @newpos\n";
#				print HAP "matches: @matches\nmaxhid: $maxhid, number of matches with this hap: $matches[$maxhid]\n" if $i == 0;
#				print HAP "hap: @{ $haps[$maxhid] }\n\n" if $i == 0 and $maxhid ne "SKIP";
				next if $maxhid eq "SKIP";	# the new snpset has no new information for any of the exisitng haps
				if ( $maxhid =~ /\d+/ ) {
					for $k ( 0 .. $#newindex ) {
						next if $haps[$maxhid][$k] eq "x";
						$haps[$maxhid][$k] = $temp[$k] if ($temp[$k] == 1 or $temp[$k] eq "x");
					}
#					print HAP "after adding snpset to this hap, hap is:\n@{ $haps[$maxhid] }\n\n" if $i == 0;
				} elsif ( $maxhid eq "NEW" ) {
					push @haps, [ @temp ];
					print "$i: Num haps: $#haps+1\n";
#					print HAP "new hap added:\n@temp\n\n" if $i == 0;
				}
			}
		}
		# check if current haps can be merged
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
# max
#
# reads in a reference for a pair of arrays
# and returns the index that has max value in first and second array (in that order).
# if max value in 2nd array = 0, return "NULL".
##################
sub max {
	my ($match, $new, $numsnps) = @_;
	@idx = ();
	$maxid = "SKIP";	# insert snpset as new hap
	$newinfo = 0;
	$conflict = 1;
	$overlap = 0;
	
	# get indices that have max match with the hap
	for my $p ( 0 .. $#{ $match } ) {
		if ($new->[$p]) {
			$newinfo += $new->[$p];
			if ($match->[$p] > 0) {		# if new information exists and matches existing  hap
				push @idx, [ $match->[$p], $new->[$p], $p ];	# conflicting haps won't have their ids stored in @idx
			}
		} else {
			$overlap = 1 if $match->[$p] == $numsnps;	# the snpset matches the exisitng hap completely
		}
		$conflict = 0 if $match->[$p] > 0;
	}

	if ( ($newinfo == 0 and $conflict == 0) or $overlap == 1 ) {
		$maxid = "SKIP";
	} elsif ($#idx > 0) {
		@sorted = sort { $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] || $a->[2] <=> $b->[2] } @idx;	# max # matches and max new sites
		# if the top two entries are same (i.e. non-unique), then confusion about which one to extend, so add a new row
#		if ($sorted[0][0] == $sorted[1][0] and $sorted[0][1] == $sorted[1][1] ) {
#			$maxid = "NEW";
#		} else {
			$maxid = $sorted[0][2];		# if no ambiguity, then add the snpset to the first one
#		}
	} elsif ( $#idx == 0 ) {
		$maxid = $idx[0][2];	# only one row
	} else {
		$maxid = "NEW";
	}
	return $maxid;
}
# End of max #
