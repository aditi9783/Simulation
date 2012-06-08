#!/usr/bin/perl -w

# reads in a hap file (output from 3maps.pl) and compares the haps.
# haps that are contained in another hap or are extendable are extended.

my $file = $ARGV[0];

my @haps;

open(FH, $file) or die "can't open file $file: $!\n";
while ( <FH> ) {
	if ($_ =~/^\d+\:\t(.+)\n/) {
		@temp = split/\s+/, $1;
		push @haps, [ @temp ];
	}
}
close FH;

my @uniqueHap;	#
my $conflict = 0;

for my $i ( 0 .. $#haps ) {
	for my $j ( 0 .. $#haps ) {
		$conflict = 0;
		for my $k ( 0 .. $#{ $haps[$i] } ) {
			my $snp = $haps[$i][$k] . $haps[$j][$k];
			if ($snp =~ /1/ and $snp =~ /x/) {
				$conflict = 1;
				last;
			}
		}
		next if $conflict == 1;
		for my $k ( 0 .. $#{ $haps[$i] } ) {
			$haps[$i][$k] = $haps[$j][$k] if ($haps[$j][$k] == 1 or $haps[$j][$k] eq "x");
			$haps[$j][$k] = 0;
			print "merging $i and $j\n";
		}
	}
}

print "Merged haps\n";
for my $i ( 0 .. $#haps ) {
	$sum = 0;
	for my $j ( 0 .. $#haps ) {
		$sum += 1 if $haps[$i][$j] == 1;
	}
	print "$i:\t@{ $haps[$i] }\n" if $sum > 0;
}
