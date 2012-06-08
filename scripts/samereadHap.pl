#!/usr/bin/perl -w

# reads in a .hapread file output by 3generatereads.pl. This file contains the read ids that contain a hap position.
# a read id with multiple hap pos would imply that the hap pos are linked. This code finds such read ids

use strict;
my ($previd, $prevline) = (-1, -1);

while ( <> ) {
	if ($_ =~ /^rid\:\s+(\d+)\,/) {
		if ($1 == $previd) {
			print $prevline;
			print $_;
			print "\n";
		}
		$previd = $1;
		$prevline = $_;
	}
}
