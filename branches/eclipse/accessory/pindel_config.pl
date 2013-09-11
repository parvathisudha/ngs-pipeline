use strict;
use Getopt::Long;
use Data::Dumper;
####### get arguments      ###
my ( $in, $bam, $id, $out, $min_insert_size );
GetOptions(
	'in=s'  => \$in,
	'bam=s' => \$bam,
	'id=s'  => \$id,
	'out=s' => \$out,
	'min_insert_size=s' => \$min_insert_size,
);

my $mean = $min_insert_size;
open IN, "<$in" or die "Can't open $in\n";
while (<IN>) {
	next unless m/^MEDIAN_INSERT_SIZE/;
	my $metrics = <IN>;
	chomp $metrics;
	my @metrics = split (/\t/, $metrics);
	$mean = $metrics[0];
	last;
}
close IN;
my $mean_rounded = sprintf "%.0f", $mean;
$mean_rounded = $min_insert_size if $mean_rounded < $min_insert_size;

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
