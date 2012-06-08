#!/usr/bin/perl -w

use strict;

# read a .reducedError file and print Precision and Recall

my (@selected, @identified, @NF);

while ( <> ) {
	if ($_ =~ /selected\:\s(\d+)\+/) {
		push @selected, $1;
	} elsif ($_ =~ /identified\:\s(\d+)/) {
		push @identified, $1;
	} elsif ($_ =~ /not\sfound\:\s(\d+)/) {
		push @NF, $1;
	}
}
print "selected:\n";
for my $p (@selected) {
	print "$p\n";
}
print "identified:\n";
for my $p (@identified) {
	print "$p\n";
}
print "haps with freq >= 0.001 that were not found:\n";
for my $p ( @NF ) {
	print "$p\n";
}
