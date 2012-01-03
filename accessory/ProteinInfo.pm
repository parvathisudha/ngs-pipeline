package ProteinInfo;
use XML::Simple;
use Data::Dumper;
use LWP::Simple;

sub new {
	my ( $class, %params ) = @_;
	my $self = {};
	my %all_params = ( %{ $params{params} }, %params );
	while ( my ( $key, $value ) = each %all_params ) {
		$self->{$key} = $value;
	}
	bless $self, $class;
	$self->{data} = $self->_download_uniprot_xml($self->{id});
	return $self;
}

sub array {
	my ( $self, $types ) = @_;
	my $data   = $self->{data};
	#print Dumper $data;
	my $result = {};
	for my $type (@$types) {
		for my $el (@$data) {
			if ( $el->{type} eq $type ) {
				$text = $el->{text}->{content} ? $el->{text}->{content} : $el->{text};
				if ( exists $result->{ $type } ) {
					$result->{ $type } =
					 $result->{ $type } . "; " . $text;
				}
				else {
					$result->{ $type } = $text;
				}
			}
		}
	}
	my @data = map{$result->{$_}} @$types;
	return \@data;
}

sub _download_uniprot_xml {
	my ( $self, ) = @_;
	return undef unless $self->{id};
	my $xml  = get 'http://www.uniprot.org/uniprot/' . $self->{id} . '.xml';
	my $data = XMLin($xml, ForceArray => [ 'gene' ],);
	my $result = $data->{entry}->{comment};
	push(@$result, { type => 'description', text => $data->{entry}->{protein}->{recommendedName}->{fullName} });
	return $result;
}

1;
