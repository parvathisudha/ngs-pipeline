use strict;
use Getopt::Long;
use Data::Dumper;
####### get arguments      ###
my ( $in, $bam, $id, $out );
GetOptions(
	'in=s'  => \$in,
	'bam=s' => \$bam,
	'id=s'  => \$id,
	'out=s' => \$out,
);

my @mean;
open IN, "<$in" or die "Can't open $in\n";
while (<IN>) {
	last unless m/^#/;
	my $mean = $1 if m/\tmean:(\d+?\.*\d+?)\t/;
	push( @mean, $mean ) if $mean;

}
close IN;

my $mean         = average( \@mean );
my $mean_rounded = sprintf "%.0f", $mean;

open OUT, ">$out" or die "Can't open $out\n";

#/data/results/Denis/S000006.20111030/S000006.20111030.bam.dedup.bam     342     S000006
print OUT join( "\t", ( $bam, $mean_rounded, $id ) ), "\n";
close OUT;

sub average {
	@_ == 1 or die('Sub usage: $average = average(\@array);');
	my ($array_ref) = @_;
	my $sum;
	my $count = scalar @$array_ref;
	foreach (@$array_ref) { $sum += $_; }
	return $sum / $count;
}
