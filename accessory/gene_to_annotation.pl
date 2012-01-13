use strict;
use Getopt::Long;
use GeneAnnotator;
use Data::Dumper;
####### get arguments      ###
my ( $in, $id_column, $uniprot, $id_type, $skip_header, $uniprot_dir);
GetOptions(
	'in=s' => \$in,
	'id_column=s'   => \$id_column,
	'uniprot=s'    => \$uniprot,
	'id_type=s' => \$id_type,
	'skip_header' => \$skip_header,
	'uniprot_dir=s' => \$uniprot_dir,
);
my $annotation = GeneAnnotator->new(id_type => $id_type, uniprot => $uniprot, uniprot_dir=> $uniprot_dir );
my $types = ['description','function', 'catalytic activity', 'tissue specificity', 'pathway', 'disease'];
open IN, $in or die "Can't open $in\n";
if(!$skip_header){
	my $header = <IN>;
	chomp $header;
	my @to_print = ($header, @$types);
	print join("\t", @to_print), "\n";		
}
while(<IN>){
	chomp;
	my @d = split /\t/;
	my $id = $d[$id_column];
	my $info = $annotation->protein_info($id, $types);
	my @to_print = ($_, @$info);
	print join("\t", @to_print), "\n";	
}
close IN;

