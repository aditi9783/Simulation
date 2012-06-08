#!/usr/bin/perl -w

use strict;

# read a .reducedError file and print Precision and Recall

my (@precision, @recall);

while ( <> ) {
	if ($_ =~ /Precision\:\s(.+)\tRecall\:\s(.+)\n/) {
		push @precision, $1;
		push @recall, $2;
	}
}
print "precision:\n";
for my $p (@precision) {
	print "$p\n";
}
print "recall:\n";
for my $p (@recall) {
	print "$p\n";
}
