package GeneAnnotator;
use Data::Dumper;

sub new {
	my ( $class, %params ) = @_;
	my $self = {};

	my %all_params = ( %{ $params{params} }, %params );
	while ( my ( $key, $value ) = each %all_params ) {
		$self->{$key} = $value;
	}
	bless $self, $class;
	if ( $self->{id_type} eq 'gene' ) {
		$self->{$self->{id_type}} = $self->read_annotation( 0, 2, $self->{gene_to_protein} );
	}
	elsif ( $self->{id_type} eq 'transcript') {
		$self->{$self->{id_type}} = $self->read_annotation( 1, 2, $self->{gene_to_protein} );
	}
	$self->{DB} = $self->read_uniprot_db;
	return $self;
}

sub get_header {
	my ( $self,) = @_;
	return $self->{DB}->{accession};
}

sub protein_info {
	my ( $self, $gene_id) = @_;
	my $protein_id = $self->gene_to_protein($gene_id);
	if($protein_id && exists $self->{DB}->{$protein_id}){
		return $self->{DB}->{$protein_id};
	}
	return undef;
}

sub gene_to_protein{
	my ( $self, $id) = @_;
	my $base = $self->{id_type};
	if(exists $self->{$base}->{$id}){
		return $self->{$base}->{$id};
	}
	else{
		return undef;
	}
}

sub read_uniprot_db {
	my ( $self,) = @_;
	my $data = {};
	my $file = $self->{uniprot_db};
	open IN, $file or die "Can't open $file\n";
	while (<IN>) {
		chomp;
		my ($id, $value) = ($1, $2) if m/(.+?)\t(.+)$/;
		next unless $id;
		$data->{ $id } = $value;
	}
	close IN;
	return $data;
}

sub read_annotation {
	my ( $self, $id_column, $an_column, $file ) = @_;
	my $data = {};
	open IN, $file or die "Can't open $file\n";
	while (<IN>) {
		chomp;
		my @d = split /\t/;
		next unless $d[$an_column];
		$data->{ $d[$id_column] } = $d[$an_column];
	}
	close IN;
	return $data;
}

1;
