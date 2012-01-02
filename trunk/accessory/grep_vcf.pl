use strict;
use Getopt::Long;
use GeneAnnotator;
####### get arguments      ###
my ( $in, $out, $regexp );
GetOptions(
	'in=s'     => \$in,
	'out=s'    => \$out,
	'regexp=s' => \$regexp,
);
open IN,  "<$in";
open OUT, ">$out";

while (<IN>) {
	if (m/^#/) {
		print OUT;
	}
	elsif (m/$regexp/) {
		print OUT;
	}
}
close IN;
close OUT;

