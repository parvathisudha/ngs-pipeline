use strict;
use Getopt::Long;
####### get arguments      ###
my ( $in, $out, $bedtools, $bed);
GetOptions(
	'in=s' => \$in,
	'out=s'   => \$out,
	'bedtools=s'    => \$bedtools,
	'bed=s' => \$bed,
);

system ("grep '#' $in > $out");
system ("$bedtools/intersectBed -u -a $in -b $bed >> $out");