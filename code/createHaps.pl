#!/usr/bin/perl -w

use strict;
use Getopt::Std;

our ($opt_e, $opt_s, $opt_h, $opt_n, $opt_m, $opt_r);
getopt( 'eshnmr' );

# set option defaults #

my $strains = 100;	# number of strains in the quasispecies
my $nhap = 1;		# number of haplotypes
my $minpos = 2;		# number of minimum positions in each haplotype
my $maxpos = 4;		# maximum number of positions in each haplotype
my $ref;		# reference seq in fasta format
my $evol = 0;		# number of sites that are different between two haplotypes under evolutionarily linked hap model

# assign user defined values
$strains = $opt_s if $opt_s;
$nhap = $opt_h if $opt_h;
$minpos = $opt_n if $opt_n;
$maxpos = $opt_m if $opt_m;
$evol = $opt_e if $opt_e;

if ($opt_r) {
	$ref = $opt_r;
} else {
	print "USAGE:\n-s\tNumber of strains in quasispecies (default = 100).";
	print "\n-h\tNumber of haplotypes (default = 1).";
	print "\n-n\tNumber of minimum positions in each haplotype (default =2).";
	print "\n-m\tNumber of maximum positions in each haplotype (default =4).";
	print "\n-r\tReference sequence in fasta format.";
	print "\n-e\tEvolve haplotypes evolutinarily: Number of sites in a new haplotype that differ from a previous hap (default = 0).\n";
	die "\n\n";
}

## Main ## 
my ($genome, $haps);	
my @nt = ("A", "T", "C", "G", "A");
my $hap;	#array of array, 1st dim: #haps, 2nd dim has hapfreq, 2D array of positions and hap bases
my $ref_freq = 0;	# frequency of reference genome
my $npos; 		# number of positions in a hap
my $prevhap;		# previous hap generated (used for creating evol linked haps)
my %seen;	# hash to store positions that are generated
my @hfreq = ();		# hap frequencies

# get the genome sequence
$genome = getSeq( $ref );

# generate hap frequencies
my $freqsum = 0;
for my $i ( 0 .. $nhap-1 ) {
	$hfreq[$i] = ( int(rand( 5 ))+1 ) / $strains;
	$freqsum += $hfreq[$i];
}
# if frequency sum > 1, normalize such that they add to 1, otherwise assign freq to reference
if ( $freqsum > 1 ) {
	for my $n ( 0 .. $#hfreq ) {
		$hfreq[$n] = sprintf( "%.4f", $hap->[$n]->[0] / $freqsum );
	} 
} else {
	$ref_freq = 1 - $freqsum;
}

# print reference genome frequency
print "\n$ref_freq";

# for each haplotype, generate hap positions
while ( $nhap > 0 ) {
	# generate hap freq, max frequency for any given haplotype is (1/20)th of #strains
	my $hapbase;		# 2D array of hap positions and hap bases

	# generate hap positions and bases 
	if ( $prevhap ) {
		$hapbase = evolHap( $prevhap ); 
		$evol = $opt_e;		# reinitialize number of sites in each new hap that differ from its parent
	} else {
		$hapbase = randomHap();
	}

	# sort $hapbase by haplotype position
	my @sorted_hapbase = sort { $a->[0] <=> $b->[0] } @{ $hapbase };
	$prevhap = \@sorted_hapbase;

	undef( $hapbase );
	$nhap--;

	# print haplotype
	print "\n$hfreq[$nhap]";
	for my $i ( 0 .. $#sorted_hapbase ) {
		print "\t$sorted_hapbase[$i][0]\t$sorted_hapbase[$i][1]";	# hap pos and base
	}
	print "\n";
}
print "\n";

################
# getSeq
#
# get reference genome sequence from
# fasta file and return it as a reference
# to an array
#
# USAGE: $seq = getSeq( $fasta_file );
################
sub getSeq {
	my $ref = shift;

	my @seq;
	open( FH, $ref ) or die "can't open file $ref:$!\n";
	while (<FH>) {
		if ($_ =~ /^[a-zA-Z]+/) {
			chomp;
			push @seq, split //;
		}
	}
	return \@seq;
}
# End of getSeq #

#################
# randomHap
#
# generate haps randomly
#
#################
sub randomHap {
	my $firsthap;
	$npos = int( rand($maxpos-$minpos+1) ) + $minpos;
	while ( $npos > 0 ) {
		my $pos = int( rand( $#{$genome}+1 ) );
		redo if $seen{$pos}++;			# redo if $pos already seen
		my $base = $genome->[$pos];
		# convert base to numeric index that matches @nt, 0:A, 1:T, 2:C, 3:G
		$base =~ tr/ATCG/0123/;
		# haplotype base is just the nt at index $base+1
		push @{ $firsthap }, [ $pos+1, $nt[$base+1] ];		# +1 to get position indexed from 1 instead of 0
#		print "base $genome->[$pos], hap pos: $pos+1, hap base: $nt[$base+1]\n";
		$npos--;
	} 
	return $firsthap;
}
# End of randomHap #

####################
# evolHap
#
# reads in previous hap and changes $evol positions in it
# to create new evolutionarily linked hap
#
####################
sub evolHap {
	my $oldhap = shift;
	my @newhap = ();
	@newhap = @{ $oldhap };

	while ( $evol > 0 ) {

		# 50/25/25 chance of creating a new hap base at this position or deleting this hap pos or creating new hap pos
		my $toss = rand();
		if ($toss <= 0.5) {	# new hap base is just the nt at index $base+1 in @nt
			my $pos = int(rand( $#newhap+1 ));		# select a hap position randomly
		#	my $base = $newhap[$pos][1];		# hap base at this position		
		 	my $base = $genome->[$pos];
			# convert base to numeric index that matches @nt, 0:A, 1:T, 2:C, 3:G
			$base =~ tr/ATCG/0123/;
			my $newbase = int( rand(4) );
			next if $newbase == $base;	# if new hap base is true base in genome, then don't do anything
#			print "changing snp base: $newhap[$pos][0], $newhap[$pos][1]. new base: $nt[$newbase]\n";
			$newhap[$pos][1] = $nt[$newbase];
		} elsif ($toss > 0.5 and $toss <= 0.75 ) { 	# add a new hap pos
			my $newpos = int(rand( $#{$genome}+1 ));
			redo if $seen{ $newpos }++;
		 	my $base = $genome->[$newpos];
			$base =~ tr/ATCG/0123/;
			my $newbase = int( rand(4) );
			redo if $newbase == $base;	# hap base should be different from true genome base
			push @newhap, [$newpos+1, $nt[$newbase]]; 
#			print "insert new snp: $newpos+1, $nt[$newbase]\n";
		} else {		# remove this hap positon or create a new one
			my $pos = int(rand( $#newhap+1 ));		# select a hap position randomly
#			print "remove snp: @{$newhap[$pos]}\n";
			splice (@newhap, $pos, 1);
		}
		$evol--;
	}
	return \@newhap;
}
# End of evolHap #
