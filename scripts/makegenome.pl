#!/usr/bin/perl -w

# generates a random genome of length l with equal probability for the four bases.

use strict;

my @base = ("A", "T", "C", "G");
my $l = $ARGV[0];

print ">genome\n";
for my $i ( 0 .. $l-1 ) {
	print "$base[ int(rand(4)) ]";
}
print "\n";
