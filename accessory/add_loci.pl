use strict;
use Getopt::Long;
use Data::Dumper;

my ($in, $out, $loci);
GetOptions( 'in=s' => \$in,
            'out=s' => \$out,
            'loci=s' => \$loci,
        );

my $loci_table = {};

open LOCI, "<$loci" or die "Can't open file: $loci\n";
while(<LOCI>){
        chomp;
        my ($chr, $start, $stop, $comment) = split /\t/;
        if(exists $loci_table->{$chr}){
        }
        else{
                $loci_table->{$chr} = [];
        }
        push(@{$loci_table->{$chr}}, [$start, $stop, $comment]);
}
close LOCI;

open IN, "<$in" or die "Can't open file: $in\n";
open OUT, ">$out"  or die "Can't open file: $out\n";
my $head = <IN>;
chomp $head;
print OUT "$head\tLOCUS\n";
while(<IN>){
        chomp;
        my @d = split /\t/;
        my $locus = get_locus($d[0], $d[1]);
        print OUT $_, "\t", $locus, "\n";
}
close OUT;
close IN;

sub get_locus{
        my ($chr, $start) = @_;
        my @cur_loci;
        for my $l (@{$loci_table->{$chr}}){
                if($start > $$l[0] && $start < $$l[1]){
                        push(@cur_loci, $$l[2]);
                }
        }
        return join(',', @cur_loci);
}
