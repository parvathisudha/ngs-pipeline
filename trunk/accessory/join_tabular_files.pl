use strict;

my ($table, $annotation) = @ARGV;

open TABLE, "<$table" or die "Can't open $table\n";
my $str = <TABLE>;
print $str;
open ANN, "<$annotation" or die "Can't open $annotation\n";

close ANN;
close TABLE;