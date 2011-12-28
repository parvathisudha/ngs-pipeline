use Data::Dumper;
use GeneAnnotator;

#my $ann = GeneAnnotator->new(transcript => 'ENST00000372726');
#$ann->get_data;
my $ann = GeneAnnotator->new(gene => 'ENSG00000124251');
$ann->get_data;