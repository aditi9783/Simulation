#!/usr/bin/perl -w

use strict;

# reads in deletions in snpsets and build a link matrix for those. Deletions in a hap would be linked together

my $file = $ARGV[0];

open( FH, $file ) or die "can't open file $file: $!\n";
my @nodes;
my @connect;

while ( <FH> ) {
	if ($_ =~ /\d+\:\s(.+)\n/) {
		my @temp = split /\s+/, $1;
		if ( $#temp == 0 ) {	# only one element, insert into @nodes is doesn't exist already
			insertNode( $temp[0] );
		}
		else {
			print "set: @temp\n";
			my $newnode = $temp[0];
			my @linked = ($newnode);
			for my $i ( 1 .. $#temp ) {
				if ($temp[$i] == $temp[$i-1]+1) {	# check if the snp indexes are continuous
					$newnode = $newnode . "-" . $temp[$i];
					print "\tnode expanded: $newnode\n";
				} else {
					insertNode( $newnode );
					print "\tnode inserted: $newnode\n";
					push @linked, $newnode;
					$newnode = $temp[$i];
				}
			}		
			insertNode( $newnode );
			print "\tnode inserted: $newnode\n";
			push @linked, $newnode;
		#	updateLinks( \@linked );
		}
	}
}
close FH;

print "Nodes are: $#nodes+1\n";
for my $n ( @nodes ) {
	print "$n\n";
}
#print "\nLinks are:\n";
#for my $i ( 0 .. $#connect ) {
#	print "$nodes[$i]: @{ $connect[$i] }\n";
#}

#########################
# insertNode
#
# inserts a node into @node if its not present already
########################
sub insertNode {
	my $newnode = shift;
	my $present = 0;
	for my $i ( 0 .. $#nodes ) {
		if ($nodes[$i] eq $newnode) {
			$present = 1;
			last;
		}
	}
	if ($present == 0) {
		push @nodes, $newnode;
	}
}
# End of insertNode #

####################
# updateLinks
#
# update the links for the nodes
####################
sub updateLinks {
	my $linkref = shift;
	for my $i ( 0 .. $#{ $linkref } ) {
		for my $j ( 0 .. $#{ $linkref } ) {
			$connect[ $nodes[ $linkref->[$i] ] ][ $nodes[ $linkref->[$j] ] ] = 1;
			$connect[ $nodes[ $linkref->[$j] ] ][ $nodes[ $linkref->[$i] ] ] = 1;
		}
	}
}
# End of updateLinks #
