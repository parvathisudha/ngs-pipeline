use strict;
use Getopt::Long;
use GeneAnnotator;
####### get arguments      ###
my ( $in, $out, $regexp, $regexp_v );
GetOptions(
	'in=s'     => \$in,
	'out=s'    => \$out,
	'regexp=s' => \$regexp,
	'regexp_v=s' => \$regexp_v,
);
open IN,  "<$in";
open OUT, ">$out";

while (<IN>) {
	if (m/^#/) {
		print OUT;
	}
	elsif (m/$regexp/) {
		unless(m/$regexp_v/){
			print OUT;
		}
	}
}
close IN;
close OUT;

