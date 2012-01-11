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
	if ( $self->{id_type} eq 'gene' ) {
		$self->{$self->{id_type}} = $self->read_annotation( 0, 2, $self->{uniprot} );
	}
	elsif ( $self->{id_type} eq 'transcript') {
		$self->{$self->{id_type}} = $self->read_annotation( 1, 2, $self->{uniprot} );
	}
	return $self;
}

sub protein_info {
	my ( $self, $gene_id, $types) = @_;
	my $uniprot_id = $self->gene_to_protein($gene_id);
	my $info = ProteinInfo->new( id => $uniprot_id);
	return $info->array($types);
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
