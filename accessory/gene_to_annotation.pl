use strict;
use Getopt::Long;
####### get arguments      ###
my ( $in, $id_column, $uniprot, $id_type);
GetOptions(
	'in=s' => \$in,
	'id_column=s'   => \$id_column,
	'uniprot=s'    => \$uniprot,
	'id_type=s' => \$id_type,
);

my $annotation = GeneAnnotator->new(gene => '1', uniprot => $uniprot,);
my $types = ['function', 'tissue specificity', 'pathway', 'disease'];
open IN, $in or die "Can't open $in\n";
my $header = <IN>;
chomp $header;
my @to_print = ($header, @$types);
print join("\t", @to_print), "\n";	
while(<IN>){
	chomp;
	my @d = split //;
	my $id = $d[$id_column];
	my $info = $annotation->protein_info($id, $types);
	my @to_print = ($annotation, @$info);
	print join("\t", @to_print), "\n";	
}
close IN;

