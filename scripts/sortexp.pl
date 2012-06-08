#!/usr/bin/perl -w

my @arr = ( [2, 3, 4], [1, 2, 5], [6, 3, 4], [6, 7, 8] );
my @sorted = sort{ $b->[0] <=> $a->[0] || $b->[1] <=> $a->[1] } @arr;

print "sorted:\n";
for my $i ( 0 .. $#sorted ) {
	print "@{ $sorted[$i] }\n";
}
