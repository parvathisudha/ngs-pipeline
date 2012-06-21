use strict;
use Getopt::Long;
use GeneAnnotator;
use Data::Dumper;
####### get arguments      ###
my ( $in, $id_column, $gene_to_protein, $id_type, $uniprot_db );
GetOptions(
	'in=s'              => \$in,
	'id_column=s'       => \$id_column,
	'gene_to_protein=s' => \$gene_to_protein,
	'id_type=s'         => \$id_type,
	'uniprot_db=s'     => \$uniprot_db,
);

my $annotation = GeneAnnotator->new(
	id_type         => $id_type,
	gene_to_protein => $gene_to_protein,
	uniprot_db      => $uniprot_db
);
my $ann_header = $annotation->get_header();
open IN, $in or die "Can't open $in\n";
my $header = <IN>;
chomp $header;
$id_column = name_to_column($header, $id_column);
warn $id_column, "\n";
print "$header\t$ann_header\n";

while (<IN>) {
	chomp;
	my @d        = split /\t/;
	my $id       = $d[$id_column];
	my $info     = $annotation->protein_info( $id);
	my @to_print = ( $_, $info );
	print join( "\t", @to_print ), "\n";
}
close IN;

sub name_to_column{
	my ($h, $name) = @_;
	my @header = split ('\t', $h);
	for(my $i = 0; $i < scalar @header; $i++){
		return $i if lc(($header[$i])) eq lc ($name);
	}
	return 0;
}