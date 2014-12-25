package Field;
use strict;
use Data::Dumper;

sub new {
        my ( $class, $id, $buffer) = @_;
        my $self = {
		id => $id,
                header2data => {},
        };
        bless $self, $class;
	$self->parse_raw_field_record($buffer);
        return $self;
}

sub parse_raw_field_record{ 
        my ( $self, $buffer) = @_;
        my @record_data = split("\n+", $buffer);
        my $data = {};
        my $field_id = "no_header";
        my $field_buffer = "";
        while(my $str = shift(@record_data)){
                if($self->is_header($str)){
			$str = make_id($str);
                        if($field_buffer){
			        $field_buffer = prepare_buffer($field_buffer);
                                $self->{header2data}->{$field_id} = $field_buffer;
                                $field_id = $str;
                                $field_buffer = "";
                        }
                        else{
                                $field_id = $str;
                        }
                }
                else{
                        $field_buffer .= $str . "\n";
                }
        }
	$field_buffer = prepare_buffer($field_buffer);
        $self->{header2data}->{$field_id} = $field_buffer;
}

sub prepare_buffer{
	my ($field_buffer) = @_;
        chomp $field_buffer;
        $field_buffer =~ s/^\s+//;
	$field_buffer =~ s/\n\s+/\n/g;
	return $field_buffer;
#       $field_buffer =~ s/(\012|\015\012?)/\\n/g;
}

sub is_header{
	my ( $self, $in ) = @_;
        my $str = $in;
        $str =~ s/\s+//g;
        $str =~ s/\///g;
	
	if($self->{id} eq 'TX'){
		return 1 if $str =~ m/^[A-Z]+$/;	
	}
	elsif($self->{id} eq 'CS'){
		return 1 if $str =~ m/^[A-Z].+:$/;
	}
	return 0;
}

sub make_id{
	my ($str) = @_;
	$str =~ s/:$//;
	$str = lc ($str);
	$str =~ s/,/ /g;	
	$str =~ s/\s+/_/g;
	$str =~ s/\//_/g;
	return $str;
}

1;

