#!/usr/bin/perl -w

# generate simulated reads for different sequencing technologies 

# 0: Illumina PE: 100 bp read length, 1% error rate, insert size 200-500 bp
# 1: Illumina MPE (mate paired end): 100 bp read length, 1% error, insert 200-8kb
# 2: Illumina SE: 100 bp read length, 1% error rate
# 3: PacBio: 3000 read length, 18% error rate

use strict;
use Getopt::Std;

our ($opt_t, $opt_n, $opt_r, $opt_h);
getopt( 'tnrh' );

## Global Var ##
# set option values
my $tech = 0;		# technology code
my $nr = 10000;		# number of reads
my $ref;

# assign user defined values
$tech = $opt_t if $opt_t;
my $hapfile = $opt_h if $opt_h;
$nr = $opt_n if $opt_n;

if ($opt_r) {
	$ref = $opt_r;
} else {
	print "USAGE:\n-r\tReference genome in fasta format\n-n\tNumber of reads to be generated (default: 10000 single/20000 paired end)\n";
	print "-t\tTechnology code:\n\t0: Illumina paired end (default)\n\t1: Illumina mate paired end\n\t2: Illumina single end\n\t3: PacBio\n";
	print "-h\tHaplotype file (optional). 1st line contains frequency of reference genome, every other line contains <hap_freq>, <hap_position1>, <hap_base1>, <hap_position2>, ..., separated by tabs. The positions are in ascending order. See tryhaps.txt for example.\n";
	die "\n\n";
}	

my @techname = ("illPE", "illMPE", "illSE", "PacBio");
my @base = ("A", "T", "C", "G");
my @rlen = (100, 100, 100, 1000);	# readlengths for diff technologies, index matches tech code
my @errpcnt = (1, 1, 1, 18);		# percent read error rate for diff tech, index matches tech code
my @rerr;				# number of sites per read that'll have substitution error
for my $i ( 0 .. $#rlen ) {
	$rerr[$i] = int( $errpcnt[$i] * $rlen[$i] / 100 ); 
}

my @ilower = (200, 200);	# lower and upper bounds for paired end data, index matched tech code
my @iupper = (500, 8000);

my ($haps, $seq, $isize, $rid, $rqual, $hposbase, $nreads, $hapflag, $errorbase);
my (@genome, @h, $baseErr, @posCov);	# baseErr is ref to 2d arrays: errors/snp at each pos in seq, posCov contains coverage at each pos
my $NHAP = 0;			# keeps track of total number of hap bases to be inserted
my ($i, $n, $pos, $b);  	# common vars to be used as index in for loops, #haps, pos in read, and hap indices

my (@ridx, @lidx, @read);	# read start and last index, identifier, actual read
				# if paired ends, then [0]-> first read, [1]-> paired read
my $debug = 0;			# prints debug statements when defined

## Main ##

getSeq();
@genome = split "", $seq;	

undef ( $seq );

my $fname = $ref . "_sim_" . $techname[$tech] . "_nr" . $nr;	
$fname = $fname . "_hfile_" . $hapfile if $hapfile;
my $fname1 = $fname . "_fir.fa";
my $fname2 = $fname . "_sec.fa" if $tech < 2;

if ( $hapfile ) {
#	open( HR, ">$fname.hapread" ) or die "can't open file to write: $!\n";	# this file contains read ids that contain hap positions
	getHaps();
}

open ( R0, ">$fname1" ) or die "can't open file to write: $!\n";
if ( $tech < 2 ) {
	open ( R1, ">$fname2" ) or die "can't open file to write: $!\n";
}

# open files to write read ids that map at each position in genome, and the errors/snp that are introduced at any given position
open (ERR, ">$fname.err") or die "can't open file to write: $!\n";
open (COV, ">$fname.coverage") or die "can't open file to write: $!\n";

# initialize base error and position read arrays
for $i ( 0 .. $#genome ) {
	push @{ $baseErr->[$i] }, 0;
}

getReads();

#print errors at each pos, and reads that map at each pos
for $i ( 0 .. $#genome ) {
	$pos = $i+1;
	print ERR "$pos:\t@{ $baseErr->[$i] }\n";
	print COV "$pos:\t$posCov[$i]\n";
}

undef (@genome);
undef ($baseErr);
undef (@posCov);

#close HR if $hapfile;
close R0;
close R1 if $tech <2;
close ERR;
close COV;

## End of main ##

###############
# getReads
#
# reads a reference genome, number of reads required
# and seq technology code (that determines read length)
# and returns a reference to a hash where key is fasta header
# for read, and key is read sequence. For paired end reads,
# data str is reference to an array of two hashes, where 
# 1st entry contains reads for one set of ends, 
# and 2nd entry contains reads for the other paired end
#
# USAGE: $reads = getReads( $ref, $num_reads, $seqtech);
###############
sub getReads {

	# generate indices for read start positions
	for $rid ( 0 .. $nr-1 ) {			# number of reads to be generated
		$ridx[0] = int(rand($#genome+1));	# read index
		$lidx[0] = $ridx[0]+$rlen[$tech]-1;
		$lidx[0] = $#genome if $#genome < $lidx[0];
		$read[0]= [ @genome[ $ridx[0] .. $lidx[0] ] ];
		for $i ( $ridx[0] .. $lidx[0] ) {
			$posCov[$i]++;
		}
		print "read $rid: idx $ridx[0]+1: @{ $read[0] }\n" if $debug;
	
		# generate paired end reads if required
		if ( $tech < 2 ) {
			# generate insert sizes for each paired end read
			$isize = int(rand($iupper[$tech] - $ilower[$tech])) + $ilower[$tech];
			$ridx[1] = $ridx[0] + $rlen[$tech] + $isize;
			$lidx[1] = $ridx[1]+$rlen[$tech]-1;	# last index for the paired read
			$lidx[1] = $#genome if $#genome < $lidx[1];
			$read[1]= [ @genome[ $ridx[1] .. $lidx[1] ] ];
			for $i ( $ridx[1] .. $lidx[1] ) {
				$posCov[$i]++;
			}
			print "paired read $rid: idx $ridx[1]+1: @{ $read[1] }\n" if $debug;
		}
		introErrors();

		print "haps left to be inserted (NHAP): $NHAP\n" if $debug;
		# insert haplotypes if requested
		if ($NHAP > 0) {
			introHaps( $rid );
		}

		# print reads
		printFasta( $rid );

		# free memory
		undef (@read);
		undef (@lidx);
		undef (@ridx);
		undef ($isize);
	}
}
# End of getReads #

###################
# getSeq
#
# reads a seq file in fasta format
# and returns the genome seq as a string
# 
# USAGE: $seq = getSeq ( $seqfile );
###################
sub getSeq {
	open(FH, $ref) or die "can't open file $ref: $!\n";
	while (<FH>) {
		chomp;
		$seq = $seq . $_ if ($_ =~ /^[a-zA-Z]+/); 
	}
	close (FH);
}
# End of getSeq #

###################
# introErrors
#
# input: reads generated by getReads
# output: reads with introduced errors
#
# USAGE: introErros( $reads, $tech );
###################
sub introErrors {
	for $i ( 0 .. $#read ) {
		for $n ( 0 .. $rerr[$tech]-1 ) {
			my $epos = int(rand( $rlen[$tech] ));	# position with error in range 0 to rlen-1 i.e. index starts at 0
			if ( $read[$i][$epos] ) {
				print "epos:$epos, orig base: $read[$i][$epos], " if $debug;
				$errorbase = $base[ int(rand(4)) ];
				redo if $errorbase eq $read[$i][$epos];				# error base should be different from true base
				$read[$i][$epos] = $errorbase;
				print "mutbase: $read[$i][$epos]\n" if $debug;
				push @{ $baseErr->[ $ridx[$i]+$epos ] }, $read[$i][$epos];	# insert the errored base at the right index
			}
		}
	}
}
# End of introErrors #

##################
# printFasta
#
# prints the reads in fasta format. Paired
# end reads are printed in two different files
# with suffixes <>_fir.fastq and <>_sec.fastq
#
# USAGE: printFastq( $reads, $tech );
###################
sub printFasta {
	$rid = shift;
	for $i ( 0 .. $#read ) {
		next if $#{ $read[$i] } == -1;  # read is of length 0, happens with paired end reads
#		$rqual = ".";			# read quality string
#		for $pos ( 0 .. $#{ $read[$i] }-1 ) {
#			$rqual = $rqual . ".";
#		}
		$pos = $ridx[$i]+1;	# read index starting from 1 (and not 0).
		if ($i == 0) {
#			print R0 "@ p0 :r $rid :idx $pos\n";
			print R0 ">p0 :r $rid :idx $pos\n";
			print R0 @{ $read[$i] };
#			print R0 "\n+\n$rqual";
			print R0 "\n";
		} elsif ($i == 1) {
#			print R1 "@ p1 :r $rid :idx $pos\n";
			print R1 ">p1 :r $rid :idx $pos\n";
			print R1 @{ $read[$i] };
#			print R1 "\n+\n$rqual\n";
			print R1 "\n";
		}
	}
}
# End of printFasta #

####################
# getHaps
#
# reads the haplotype file and returns the hap info as 
# array of array where for each hap (row) contains a 
# 2d array where each row = [hap pos, hap base, #reads with this hap]
#
# USAGE: $haps = getHaps( $hapfile );
####################
sub getHaps {
	my $coverage = $nr * $rlen[$tech] / $#genome;
	$coverage = $coverage * 2 if $tech < 2;		# number of reads are doubled for paired end methods
	# store haps in $haps, 1st value is the hap frequency
	open( HAP, $hapfile ) or die "can't open file $hapfile: $!\n";
	while ( <HAP> ) {
		if ($_ =~ /.+\t\d+/) {
			chomp;
			@h = split /\t/;	# haplotype in each line of hapfile
			$nreads = int( $h[0] * $coverage )+1;	# numer of reads that should have this haplotype, at least once
#			print HR "haplotype: @h\t#reads for each pos: $nreads\n";
			for ( $i=1; $i<=$#h; $i +=2 ) {
				push @{ $hposbase }, [$h[$i]-1, $h[$i+1], $nreads];	# hap pos (-1 to start index from 0), base, and #reads it should be in
				print "\tpos $h[$i], base $h[$i+1]" if $debug;
				$NHAP += $nreads;
				print "\tNHAP $NHAP\n" if $debug;
				last if $i >= $#h;
			}
			push @{ $haps }, $hposbase;
		}
		undef (@h);
		undef ($hposbase);
		undef ($nreads);
	}
}
# End of getHaps #

####################
# introHaps
#
# reads in a file containing haplotypes and their
# respective frequencies and intoduces those 
# haplotypes in the reads
#
# USAGE: introHaps( $reads, $hapfile, $index, $tech, $seqfile );
####################
sub introHaps {
	$rid = shift;
	for $n ( 0 .. $#{ $haps } ) {
		$hapflag = 0;
		# check if this read covers a hap. A read (or paired read) is allowed to have one haplotype only, but can contain all pos of this hap
		for $i ( 0 .. $#read ) {		# for each read
			for $b ( 0 .. $#{ $haps->[$n] } ) { 	# for all bases in this hap
				next if $haps->[$n]->[$b]->[2] == 0;	# skip if this hap base has been inserted in required number of reads
				
				# if hap position lies in read, insert hap base
				if ( $haps->[$n]->[$b]->[0] >= $ridx[$i] and $haps->[$n]->[$b]->[0] <= $lidx[$i] ) {
					$read[$i]->[ $haps->[$n]->[$b]->[0] - $ridx[$i] ] = $haps->[$n]->[$b]->[1];
					$haps->[$n]->[$b]->[2]--;	# decrease counter for # reads for this hap pos
					$NHAP--;			# decrease overall # hap counts
					$hapflag++;
					print "\t\thap found: rid $rid, happos $haps->[$n]->[$b]->[0], NHAP $NHAP\n" if $debug;
					$pos = $haps->[$n]->[$b]->[0] + 1;
#					print HR "rid: $rid, hap pos: $pos, base: $haps->[$n]->[$b]->[1]\n";
					push @{ $baseErr->[ $haps->[$n]->[$b]->[0] ] }, $haps->[$n]->[$b]->[1];
				}
			}
		}
		last if $hapflag > 0;	# allow only one haplotype for each read
	}	
}
# End of introHaps #
