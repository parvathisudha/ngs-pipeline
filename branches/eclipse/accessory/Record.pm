package Record;
use strict;
use Data::Dumper;
use Field;

sub new {
        my ( $class, $buffer) = @_;
        my $self = {
		id => undef,
                field2data => undef,
        };
        bless $self, $class;
	$self->parse_raw_omim_record($buffer);
        return $self;
}

sub get_id{
	my ($self) = @_;
	return $self->{id};
}

#*RECORD*
#*FIELD* NO
#100050
#*FIELD* TI
#100050 AARSKOG SYNDROME, AUTOSOMAL DOMINANT
#*FIELD* TX
#
#DESCRIPTION
#
#Aarskog syndrome is characterized by short stature and facial, limb, and
#genital anomalies. One form of the disorder is X-linked (see 305400),
#but there is also evidence for autosomal dominant and autosomal
#recessive (227330) inheritance (summary by Grier et al., 1983).
#
#CLINICAL FEATURES
#
#Grier et al. (1983) reported father and 2 sons with typical Aarskog
#syndrome, including short stature, hypertelorism, and shawl scrotum.
#Stretchable skin was present in these patients.

sub parse_raw_omim_record{ 
        my ( $self, $buffer) = @_;
	my @record_data = split("\n+", $buffer);
	my $data = {};	
	shift(@record_data);
	$self->{id} = shift(@record_data);

	my $field_id = "";
	my $field_buffer = "";
        while(my $str = shift(@record_data)){
                if($str =~ m/\*FIELD\*.(..)/){
                        if($field_buffer){
                                my $field = Field->new($field_id, $field_buffer);
				$self->{field2data}->{$field_id} = $field;
                                $field_id = $1;
                                $field_buffer = "";
                        }
			else{
				$field_id = $1;
			}
                }
                else{
                        $field_buffer .= $str . "\n";
                }
        }
	my $last_field = Field->new($field_id, $field_buffer);
	$self->{field2data}->{$field_id} = $last_field;
}
sub get_field{
	my ( $self,$query ) = @_;
	my ( $field,$header ) = split (",", $query);
	my $data = $self->{field2data}->{$field}->{header2data}->{$header};
	return $data;
}

1;

