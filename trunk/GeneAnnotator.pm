package GeneAnnotator;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Registry;
use Data::Dumper;

sub new {
	my ( $class, %params ) = @_;
	my $self = {};

	my %all_params = ( %{ $params{params} }, %params );
	while ( my ( $key, $value ) = each %all_params ) {
		$self->{$key} = $value;
	}
	bless $self, $class;
	my $registry = 'Bio::EnsEMBL::Registry';

	$registry->load_registry_from_db(
		-host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
		-user => 'anonymous'
	);
	$self->{registry} = $registry;
	my $gene_adaptor = $registry->get_adaptor( 'Human', 'Core', 'Gene' );
	$self->{gene_adaptor} = $gene_adaptor;
#	my $transcript_adaptor =
#	  $registry->get_adaptor( 'Human', 'Core', 'Transcript' );
#	$self->{transcript_adaptor} = $transcript_adaptor;

	return $self;
}

sub get_data {
	my ($self) = @_;
	#my $transcript = $self->{transcript_adaptor}->fetch_by_stable_id($self->{transcript});
	my $gene = $self->{gene_adaptor}->fetch_by_stable_id($self->{gene});
	print $gene->display_id;
	print Dumper $gene->get_all_DBLinks() ;
}

1;
