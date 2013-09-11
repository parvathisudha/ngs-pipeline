use strict;
use Data::Dumper;
my $coverage = {};
my $length = {};
open IN, $ARGV[0] or die "Can't open input file\n";
while(<IN>){
        chomp;
        next if m/Y\t/;
        next if m/X\t/;
        next if m/genome\t/;
#       chrM    57      1       16571   6.03464e-05
        my @data = split("\t", $_);
        $length->{$data[0]} = $data[3];
        update_coverage(@data[1,2]);
}
close IN;

my $genome_length = get_genome_length();

my $seq_num = 0;

for (sort { $a <=> $b } keys %$coverage){
        $seq_num += $coverage->{$_};
        print $_, "\t", $seq_num/$genome_length, "\t", $seq_num, "\t", $genome_length, "\n";
}

#print "GENOME LENGTH: $genome_length\n";

sub get_genome_length{
        my $num = 0;
        for(keys %$length){
                $num += $length->{$_};
        }
        return $num;
}

sub update_coverage{
        my ($cov, $num) = @_;
        if(exists $coverage->{$cov}){
                $coverage->{$cov} = $coverage->{$cov} + $num;
        }
        else{
                $coverage->{$cov} = $num;
        }
}