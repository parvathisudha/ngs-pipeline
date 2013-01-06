use strict;
use Getopt::Long;
####### get arguments      ###
my ( $in, $out, $regexp, $regexp_v );
GetOptions(
	'in=s'     => \$in,
	'out=s'    => \$out,
	'regexp=s' => \$regexp,
	'regexp_v=s' => \$regexp_v,
);
open IN,  "<$in" or die "Can't open $in file!\n";
open OUT, ">$out" or die "Can't open $out file!\n";

my $head = <IN>;
print OUT $head;

while (<IN>) {
	if (m/$regexp/) {
		if($regexp_v){
			unless(m/$regexp_v/){
				print OUT;
			}			
		}
		else{
			print OUT;
		}
	}
}
close IN;
close OUT;

