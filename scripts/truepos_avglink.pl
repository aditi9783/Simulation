#!/usr/bin/perl -w

# reads in .avglink file and output of redErr.pl that has unique hap pos and outputs the avg links for the positions that overlap in these two files

use strict;

my $errout = $ARGV[0];
my $avglink = $ARGV[1];

my @unique_hpos = ();
my (@uhpos_alink, @other_alink);

open (ERR, $errout) or die "can't open file $errout: $!\n";
open (AVG, $avglink) or die "can't open file $avglink: $!\n";

while ( <ERR> ) {
	if ($_ =~ /^(\d+)\n/) {
		push @unique_hpos, $1;
	}
}
close ERR;
my $flag = 0;

while ( <AVG> ) {
	$flag = 0;
	if ($_ =~ /^(\d+)\:\s+(.+)\n/) {
		for my $pos ( @unique_hpos ) {
			if ($1 == $pos) {
				push @uhpos_alink, [$1, $2];
				$flag = 1;
			}
		}
		push @other_alink, [$1, $2] if $flag == 0;
	}	
}

print "Unique hap pos avglinks:\n";
for my $i ( 0 .. $#uhpos_alink ) {
	print "@{ $uhpos_alink[$i] }\n";
}

print "Other pos avglinks:\n";
for my $i ( 0 .. $#other_alink ) {
	print "@{ $other_alink[$i] }\n";
}
