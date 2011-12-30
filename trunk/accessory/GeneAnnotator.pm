package GeneAnnotator;
use ProteinInfo;
use Data::Dumper;

sub new {
	my ( $class, %params ) = @_;
	my $self = {};

	my %all_params = ( %{ $params{params} }, %params );
	while ( my ( $key, $value ) = each %all_params ) {
		$self->{$key} = $value;
	}
	bless $self, $class;

	if ( $self->{gene} ) {
		my $gene2uniprot = read_annotation( 0, 2, $ensemble_to_uniprot );
		$self->{gene2uniprot} = $gene2uniprot;
	}
	if ( $self->{transcript} ) {
		my $transcript2uniprot = read_annotation( 1, 2, $ensemble_to_uniprot );
		$self->{transcript2uniprot} = $transcript2uniprot;
	}
	return $self;
}

sub protein_info {
	my ( $self, $gene_id, $types) = @_;
	my $uniprot_id = $self->gene_to_protein($gene_id);
	my $info = ProteinInfo->new( id => $uniprot_id);
	return $info;
}

sub gene_to_protein{
	my ( $self, $id) = @_;
	if(exists $self->{gene2uniprot}->{$id}){
		return $self->{gene2uniprot}->{$id};
	}
	else{
		return undef;
	}
}

sub read_annotation {
	my ( $self, $id_column, $an_column, $file ) = @_;
	my $data = {};
	open IN, or die "Can't open $file\n";
	while (<IN>) {
		chomp;
		my @d = split //;
		next unless $d[$an_column];
		$data->{ $d[$id_column] } = $d[$an_column];
	}
	close IN;
	return $data;
}

1;
