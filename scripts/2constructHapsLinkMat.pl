#!/usr/bin/perl -w

# constructs haplotypes from the link matrix output by mapSNPs.pl (saved in .links file)

use strict;

my $lfile = $ARGV[0];
my $ethres = 0.2;		# threshold avg # of links: if value < this threshold, likely an error postion (error positions will have lower avg # links than true haps)

my (@errpos, @happos, @repeatpos, @links, @index);	# error positions, positions that have more than one occurence (with a different base), link matrix, snp at a given index
my (@min, @avglinks, @nonzeroelt);	# min number of links per row, freq dist for @min, avg number of links for each row, number of non-zero elements at each row. 
my (@haps, @sorted_min, @postrack);
my $prevpos = 0;
my $idx = 0;
my ($i, $j, $m, $eidx, $minval, $minidx, $minnonzero, @temp, $flag, $sum);					# temporary var

##### MAIN ########
readLinkFile();			# reads link matrix file and populates @min, @nonzeroelt, @errpos, @happos, and @links
#removeErr();			# remove error position from link matrix and adjust the link counts accordingly
constructHaps();

# ENd of MAIN #####

####################
# readLinkFile 
#
# reads the link file
####################
sub readLinkFile {
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
#	print "81st row: @{ $links[81] }\n";
#	updateMinAndNonzeroelt();
#	printProgress();	# prints @max, @nonzeroelt, @links	
	print "num hap pos: $#happos\n";
	for $i ( 0 .. $#happos ) {
		print "@{ $index[ $happos[$i] ] }\n";
	}
}

# End of readLinkFile #

####################
# removeErr
#
# removes the rows and cols in link matrix that
# correspond to the error positions. Also removes
# the link matrix counts for those positions
####################
sub removeErr {
	my $errors = errsum();
	print "Numbr of errors: $errors\n";
	while ( $errors ) {
		getLinkages( \@errpos );
		updateMinAndNonzeroelt();
		$errors = errsum();
		print "Number of errors left: $errors\n";
		printProgress() if $errors < 100;
	}
#	for $eidx ( @errpos ) {
#		# set remaining values in error positions to zero
#		for $i ( 0 .. $#{ $links[$eidx] } ) {
#			$links[$eidx][$i] = 0;
#			$links[$i][$eidx] = 0;
#		}
#	}
}
# removeErr #

#####################
# errsum
#
# sum of links at error positions.
#####################
sub errsum {
	my $sum = 0;
	for $i ( @errpos ) {
		for $j ( 0 .. $#{ $links[$i] } ) {
			$sum += $links[$i][$j];
		}
	}
#	print "sum of error positions: $sum\n";
	return $sum;
}
# ENd of errsum #

#####################
# getNumHaps
#
# get number of haps from the link matrix
#####################
sub getNumHaps {
	my @sorted_nze = sort {$b <=> $a} @nonzeroelt;
	my ($nd, @newpos, $nhap);
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
		$nd = distinctValArray( \@newpos );		# get number of disict values among the new positions
		$nhap += $nd;

		$flag = 0;
		for $j ( 0 .. $#postrack ) {
			$flag = 1 if $postrack[$j] == 0;	# at least one position has not been accounted for
		}
		last if $flag == 0;				# no more new positions to track
	}
	return $nhap;
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
# getLinkages
#
# reads in link matrix and constructs haps
#####################
sub getLinkages {
	my ($idxarray, $hapflag) = @_;
	($minidx, $minval) = selectIndex( $idxarray );		# minidx: select smallest index that has min link value as well as smallest number of non zero elements
							# minval: the least number of links for the positions that make this hap. This is an estimate of hap freq
#	die "no more min index to be found\n" if $minidx eq "X";
	my @linkageidx = ( $minidx );
#	print "inside getLinkages: minidx: $minidx, its links: @{ $links[$minidx] }\n";
	# the indices that this index links to forms a haplotype, and the min link value in the row is an estimate of its abundance
	for $i ( 0 .. $#{ $links[$minidx] } ) {
		if ($links[$minidx][$i] > 0) {
			push @linkageidx, $i;
		}
	}
#	print "hap ids: @linkageidx\n";
	push @haps, [$minidx, $minval, \@linkageidx] if $hapflag;	# the function is called to construct haps and not to remove errors
#
	my @row1 = @{ $links[$minidx] };

	# update link matrix to remove the $minval links from concerned indices
	for $i ( 0 .. $#linkageidx ) {
		for $j ( $i+1 .. $#linkageidx ) {
#			print "remove $minval from $linkageidx[$i],$linkageidx[$j]\n";
			my $val = $links[ $linkageidx[$i] ][ $linkageidx[$j] ]; 
			$links[ $linkageidx[$i] ][ $linkageidx[$j] ] -= $minval;
			$links[ $linkageidx[$j] ][ $linkageidx[$i] ] = $links[ $linkageidx[$i] ][ $linkageidx[$j] ];
			if ($links[ $linkageidx[$i] ][ $linkageidx[$j] ] < 0) {
				print "minidx: $minidx, minval: $minval, index pair: $linkageidx[$i], $linkageidx[$j], pair val: $val\nlinkages: @linkageidx\n";
				print "minidx row:\n@row1\n";
				die;
				
			}
		}
	}
	updateMinAndNonzeroelt();
#	printProgress();	# prints @max, @nonzeroelt, @links	
}
# End of getLinkages #

######################
# constructHaps
#
# reads in indices for non error positions and constructs haps
######################
sub constructHaps {
#	my $hapidx = shift;
	my $nhap = getNumHaps();			# computes number of haplotypes from the link matrix
	
	print "Number of haps is: $nhap\n";
	while ( $nhap > 0 ) {
		getLinkages( \@happos, 1 );		# Each call identifies one hap and stores that into global array @haps
		$nhap--;
		if ($nhap == 0) {
			printHaps();
			die "\n";
		}	
	}
}
# End of constructHaps #

#####################
# selectIndex
#
# select smalled index with min link value and also least number of non zero elts
#####################
sub selectIndex {
	my $idxarray = shift;
	my $minlink = min( $idxarray, \@min );
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
	my ($idxarray, $aref) = @_;
	$m = 10000;
	for $i ( @{ $idxarray } ) {
		next if $aref->[$i] == 0;	# looking for a nonzero min value
		$m = $aref->[$i] if $aref->[$i] < $m; 
	}
#	print "inside min , index array: @{ $idxarray }, aray is @{ $aref }, min val is $m\n";
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
#	print "81, 155: $links[81][155]\t81, 480: $links[81][480]\t155, 480: $links[155][480]\n";
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
		print "index with min freq = $haps[$i][0], freq: $haps[$i][1]: @{ $haps[$i][2] }\n";
	}
}
# End of printHaps #
