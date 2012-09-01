use strict;
use Data::Dumper;
use Getopt::Long;
####### get arguments      ###
my ( $input_file, $result_file, $normalization, );
GetOptions(
	'input=s'       => \$input_file,
	'result=s'       => \$result_file,
	'normalization=s' => \$normalization,
);

my $raw_distribution = get_distribution_from_file($input_file);
my $norm_distribution = get_distribution_from_file($normalization);
my $corrected_distribution = correct_distribution($raw_distribution, $norm_distribution);
open DISTRIBUTION, ">$result_file" or die "Can't open $result_file\n";
print_distribution($corrected_distribution, \*DISTRIBUTION);
close DISTRIBUTION;

sub get_distribution_from_file{
	my ($file) = @_;
	my $d = {};
	open IN, "<$file" or die "Can't open $file\n";
	my $h = <IN>;
	while(<IN>){
		chomp;
		my @a = split /\t/;
		$d->{$a[0]} = $a[1];
	}
	close IN;
	return $d;
}
sub correct_distribution{
        my ($d, $n) = @_;
	my $cutoff = 13;
	die "Different number of bins in distributions!\n" unless scalar(keys %$d) == scalar(keys %$n);
	my $ratios_sum = 0;
	for my $key(sort { $a <=> $b } keys %$n){
		my $ratio = $n->{$key} / $d->{$key};
		$ratios_sum += $ratio;
		last if $key == $cutoff;
	}
	my $norm_coeff = $ratios_sum / ($cutoff + 1);
        for my $key(sort keys %$d){
                $d->{$key} = $norm_coeff * $d->{$key};
        }
	return $d;
}
sub print_distribution{
	my ($d, $OUT) = @_;
	print $OUT "REPEAT_LENGTH\tREADS\n";
	for ( sort { $a <=> $b } keys %$d ) {
		print $OUT $_, "\t", $d->{$_}, "\n";
	}
}
