#!/usr/bin/perl -w

# constructs haplotypes from the link matrix output by mapSNPs.pl (saved in .links file)

use strict;

my $lfile = $ARGV[0];
my $ethres = 0.2;		# threshold avg # of links: if value < this threshold, likely an error postion (error positions will have lower avg # links than true haps)

my (@errpos, @happos, @repeatpos, @links, @index);	# error positions, positions that have more than one occurence (with a different base), link matrix, snp at a given index
my (@min, %minhist, @avglinks, @nonzeroelt);	# min number of links per row, freq dist for @min, avg number of links for each row, number of non-zero elements at each row. 
my (@haps, @sorted_min, $nhap, @postrack);
my $prevpos = 0;
my $idx = 0;
my ($i, $j, $m, $eidx, $minval, $minidx, $minnonzero, @temp, $flag, $sum);					# temporary var

open( FH, $lfile ) or die "can't open file $lfile: $!\n";

while (<FH>) {
	if ($_ =~ /^(\d+)\:\s(\d+)\s([A-Z])\n/) {
		$index[$1] = [$2, $3];		# $1 is index, $2 is snp pos, and $3 is snp base
		push @repeatpos, ($1-1, $1) if $prevpos == $2;	# store indices for repeated positions, as they likely are not error positions
		$prevpos = $2;
	} elsif ($_ =~ /^\d+\s\d+/) {
		chomp;
		@temp = split /\s+/, $_;

		$sum = 0;	
		for $i ( @temp ) {
			$sum += $i;
		}
		$avglinks[$idx] = $sum/($#temp+1);

		$flag = 0;
		if ( $avglinks[$idx] <= $ethres ) {	# error position?
			for $i ( 0 .. $#repeatpos ) {
				$flag = 1 if $idx == $repeatpos[$i];	# flag = 1 means this is a repeated position and hence likely not an error
			}
			if ( $flag == 0 ) {	# error position
				push @errpos, $idx;
				$postrack[$idx] = 1;	# error position tracked
			} else {
				push @happos, $idx;	# true hap pos index
			}
		} else {
			push @happos, $idx;
		}

		@{ $links[$idx] } = @temp;
		$idx++;	
	}
}

removeErr();			# remove error position from link matrix and adjust the link counts accordingly
updateMinAndNonzeroelt();	# update @min and @nonzeroelt
getNumHaps();			# computes number of haplotypes from the link matrix
printProgress();		# prints @min, @nozeroelt, @links
print "Number of haps is: $nhap\n";
constructHaps();

####################
# removeErr
#
# removes the rows and cols in link matrix that
# correspond to the error positions. Also removes
# the link matrix counts for those positions
####################
sub removeErr {
	for $eidx ( @errpos ) {
#		my (@row1, @row2);
		my @nze = ();		# indices of non zero elts in this row
		for $i ( 0 .. $#{ $links[$eidx] } ) {
			push @nze, $i if $links[$eidx][$i] > 0;
		}
		print "nze: @nze\nerror row:\n@{ $links[$eidx] }\n";
		for $i ( 0 .. $#nze ) {
			$minval = $links[ $eidx ][ $nze[$i] ];
			next if $minval == 0;
			for $j ( $i+1 .. $#nze ) {
				$minval = $links[ $eidx ][ $nze[$j] ] if $minval > $links[$eidx][ $nze[$j] ];
				next if $minval == 0;
				$links[$nze[$i]][$nze[$j]] -= $minval;
				$links[$nze[$j]][$nze[$i]] = $links[$nze[$i]][$nze[$j]];		# to keep matrix symmetric
				if ($links[ $nze[$i] ][ $nze[$j] ] < 0) {
					print "error row:\n@{ $links[$eidx] }\n";
					die "error index: $eidx, index pair: $nze[$i], $nze[$j],  minval: $minval\nrow $nze[$i]:\n@{ $links[$nze[$i]] }\nrow $nze[$j]:\n@{ $links[$nze[$j]] }\n";
				}
				$links[$eidx][$nze[$i]] -= $minval;
				$links[$nze[$i]][$eidx] = $links[$eidx][$nze[$i]];
				$links[$eidx][$nze[$j]] -= $minval;
				$links[$nze[$j]][$eidx] = $links[$eidx][$nze[$j]];
			}
		}
		# set remaining positions to zero
		for $i ( 0 .. $#{ $links[$eidx] } ) {
			$links[$eidx][$i] = 0;
			$links[$i][$eidx] = 0;
		}
#		print "Error positin index: $eidx\n";
#		@temp = @{ $links[$eidx] };	# number of links of this error position index with other indices in link matrix
#		for $i ( 0 .. $#temp ) {        # $i is the index of snp that share read with this error pos
#			next if $temp[$i] == 0;
#			$minval = $temp[$i];
#			for $j ( $i+1 .. $#temp ) {
#				next if $temp[$j] == 0;
#				$minval = $temp[$j] if $temp[$j] < $minval;
#				@row1 = @{ $links[$i] };
#				@row2 = @{ $links[$j] };
				# do error reduction only for non error positions, error positions will be marked with zeros anyways
#				if ($avglinks[$i] > $ethres and $avglinks[$j] > $ethres ) {
#					$links[$i][$j] -= $minval;
#					$links[$j][$i] = $links[$i][$j];		# to keep matrix symmetric
#				}
				
#				if ($links[$i][$j] < 0) {
#					print "row i $i:\n@row1\njth elt: $row1[$j]\n";
#					print "row j $j:\n@row2\nith elt: $row2[$i]\n";	
#					die "error index: $eidx, index pair: $i, $j,  minval: $minval\nlink row:\n@temp\n";
#				}
#				@row1 = ();
#				@row2 = ();
#				print "Error index pair: $i, $j, minval: $minval\n";

				# set error position row and column to contain zeros
#				$links[$eidx][$i] = 0;
#				$links[$i][$eidx] = 0;
#				$links[$eidx][$j] = 0;
#				$links[$j][$eidx] = 0;
#			}
#		}
	}
}
# removeErr #

#####################
# getNumHaps
#
# get number of haps from the link matrix
#####################
sub getNumHaps {
	my @sorted_nze = sort {$b <=> $a} @nonzeroelt;
	my ($nd, @newpos);
	for $i ( 0 .. $#sorted_nze ) {
		# get index for row that has next max number of non zero elements
		for $j ( 0 .. $#nonzeroelt ) {
			$idx = $j if $nonzeroelt[$j] == $sorted_nze[$i];
		}
		@temp = @{ $links[$idx] };
		# get link values at positions that have not been visited before. Out of these link values, get # of distinct values (=number of new haps that this row is part of)
		for $j ( 0 .. $#temp ) {
			next if $postrack[$j];		# postrack[$j] is defined: this position has already been accounted for
			push @newpos, $temp[$j];
			$postrack[$j] = 1;
		}
		$nd = distinctValArray( \@newpos );	# get number of disict values among the new positions
		$nhap += $nd;

		$flag = min( @postrack );
		last if $flag == 1;			# no more new positions to track
	}
}
# End of getNumHaps #

#####################
# distinctValArray
#
# reads in a list and returns number of disintct values (other than zero) in it
#####################
sub distinctValArray {
	my $aref = shift;
	my @sorted = sort{ $a<=> $b } @$aref;
	my $prevval = 0;
	my $num = 0;
	for my $k ( 0 .. $#sorted ) {
		next if $sorted[$k] == 0;
		$num++ if ($sorted[$k] != $prevval);
		$prevval = $sorted[$k];
	}
	return $num;
}
# End of distinctValArray #

#####################
# constructHaps
#
# reads in link matrix and constructs haps
#####################
sub constructHaps {
	($minidx, $minval) = selectIndex();		# minidx: select smallest index that has min link value as well as smallest number of non zero elements
							# minval: the least number of links for the positions that make this hap. This is an estimate of hap freq
	my @hapidx = ( $minidx );
#	print "inside constructHaps: minidx: $minidx, its links: @{ $links[$minidx] }\n";
	# the indices that this index links to forms a haplotype, and the min link value in the row is an estimate of its abundance
	for $i ( 0 .. $#{ $links[$minidx] } ) {
		if ($links[$minidx][$i] > 0) {
			push @hapidx, $i;
		}
	}
#	print "hap ids: @hapidx\n";
	push @haps, [$minval, \@hapidx];

	# update link matrix to remove the $minval links from concerned indices
	for $i ( 0 .. $#hapidx ) {
		for $j ( $i+1 .. $#hapidx ) {
#			print "remove $minval from $hapidx[$i],$hapidx[$j]\n";
			$links[ $hapidx[$i] ][ $hapidx[$j] ] -= $minval;
			$links[ $hapidx[$j] ][ $hapidx[$i] ] = $links[ $hapidx[$i] ][ $hapidx[$j] ];
		}
	}

	updateMinAndNonzeroelt();
	$nhap--;
#	printProgress();	# prints @max, @nonzeroelt, @links	
	if ($nhap == 0) {
		printHaps();
		die "\n";
	}	
	constructHaps();	# recursive call
}
# End of constructHaps #

#####################
# selectIndex
#
# select smalled index with min link value and also least number of non zero elts
#####################
sub selectIndex {
	my $minlink = min( @min );
	$m = "X";
	$minval = 100000;
	# get indices for all the snps that have this min val
	for $idx ( 0 .. $#min ) {
		next if $min[$idx] != $minlink;
		next if $nonzeroelt[ $idx ] == 0;		# all elements in this row are zero
		if ($nonzeroelt[ $idx ] < $minval) {
			$m = $idx;
			$minval = $nonzeroelt[ $idx ];
		}
#		print "inside selectIndex: idx: $idx, num of non zero elt: $minval\n";
#		print "inside selectIndex: min value: $minlink, index returned: $m\n";
	}
	return ($m, $minlink);
}
# End of selectIndex #

#####################
# updateMinAndNonzeroelt
#
# updates arrays that store min values of each row
# and number of nonzero elt in each row based on the 
# current link matrix
#####################
sub updateMinAndNonzeroelt {
	#update @min: store min link values for each row, and @nonzeroelt: number of nonzero elts in the row
	for $i ( 0 .. $#links ) {
		$nonzeroelt[$i] = 0;
		$min[$i] = 100000;
		for $j ( 0 .. $#{ $links[$i] } ) {
			if ($links[$i][$j] > 0) {
				$min[$i] = $links[$i][$j] if $links[$i][$j] < $min[$i];
				$nonzeroelt[$i]++;
			}
		}
	}
}
# End of updateMinAndNonzeroelt #

#####################
# min
#
# returns minimum value from a list
#####################
sub min {
	my @val = @_;
	my $m = 10000;
	for $i ( 0 .. $#val ) {
		next if $val[$i] == 0;	# looking for a nonzero min value
		$m = $val[$i] if $val[$i] < $m; 
	}
	return $m;
}
# End of min #

####################
# printProgress
#
# prints @max, @nozeroelt, and @links
####################
sub printProgress {
	print "Link matrix:\n";
	for $i ( 0 .. $#links ) {
		print "$i: @{ $index[$i] }: @{ $links[$i] }\n";
	}
	print "\nMin value and number of non zero elts:\n";
	for $i ( 0 .. $#min ) {
		print "$i: $min[$i]\t$nonzeroelt[$i]\n";
	}
}
# End of printProgress #

#####################
# printHaps
#
# reads in 2d array with hap indices and hap abundance
#####################
sub printHaps {
	print "Haplotypes are:\n";
	for $i ( 0 .. $#haps ) {
		print "$haps[$i][0]: @{ $haps[$i][1] }\n";
	}
}
# End of printHaps #
