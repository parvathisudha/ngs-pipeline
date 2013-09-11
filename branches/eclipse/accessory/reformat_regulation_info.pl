use strict;
use Getopt::Long;
use Data::Dumper;
####### get arguments      ###
my ( $in, $skip_header, $eff_column);
GetOptions(
	'in=s' => \$in,
	'eff_column=s' => \$eff_column,
	'skip_header' => \$skip_header,
);
open IN, $in or die "Can't open $in\n";
if(!$skip_header){
	my $header = <IN>;
	print $header;
}

while(<IN>){
	chomp;
	my @d = split /\t/;
	
	my $eff = $d[$eff_column];
	my @eff = split (/\),/ , $eff);
	my @regs;
	for my $e(@eff){
		push (@regs, $1) if $e =~ m/REGULATION\[(.+?)\]/;
	}
	$d[$eff_column] = join(",",@regs);
	print join("\t",@d), "\n";
}
close IN;

