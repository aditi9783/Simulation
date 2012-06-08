#!/usr/bin/perl -w

# compares haps that are predicted vs true haps from the files that contain those

use strict;

my $truehapfile = $ARGV[0];
my $predhapfile = $ARGV[1];
my $genomefile = $ARGV[2];

my $genome = readSeq( $genomefile );
my ($truehap, $truehapf) = readHapfile( $truehapfile );
my ($predhap, $predhapf) = readHapfile( $predhapfile );



#########################
# readHapfile
#
# reads hap file and returns hap freq and hap base-positions 
# as references to two arrays
#########################
sub readHapfile {
	my $filename = shift;
	my ($hap, $hapfreq);

	open (FH, $filename) or die "can't open file $filename : $!\n";
	while ( <FH> ) {
		next if $_ =~ /\d+.{0,10}\n/;
		chomp;
		my @h = split /\t/;
		push @{$hapfreq}, $h[0];
		my @temp = ();
		for (my $i = 1; $i < $#h; $i = $i+2) {
			if ($genome->[ $h[$i]-1 ] ne $h[$i+1]) {
				push @temp, $h[$i] . "-" . $h[$i+1];
			} else {
				print "snp matches true genome base: $h[$i] $h[$i+1]\n";
			}
		}
		push @{ $hap }, [ @temp ];
	}
	close FH;
	return ($hap, $hapfreq);
}
# End of readHapfile #

######################
# readSeq
#
# reads a seqfile (fasta format) and returns the genome 
# sequence as a ref to array
######################
sub readSeq {
	my $seqfile = shift;
	my $seq;
	my @seqarr;
	open(FH, $seqfile) or die "can't open sequence file: $seqfile: $!\n";
	while (<FH> ) {
		next if $_ =~ /\>/;
		chomp;
		$seq = $seq . $_ if $_ =~ /[A-Z]+/i;
	}
	@seqarr = split //, $seq;
	return \@seqarr;
}
# End of readSeq

