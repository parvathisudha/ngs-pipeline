package Lane;
use strict;

sub new {
	my ( $class, $dir, $num, $archived, $paired ) = @_;

	my $self = {
		DIR            => $dir,
		NUM            => $num,
		ARCHIVED       => $archived,
		PAIRED       => $paired,
	};
	bless $self, $class;
	return $self;
}
#
#sub forward_align_id {
#	my $self = shift;
#	return _get_id($self->get_read_id(1));
#}

#sub _get_id{
#	my $file = shift;
#	$file =~ s/\//_/g;
#	$file =~ s/\://g;
#	return $file;
#}

#sub reverse_align_id {
#	my $self = shift;
#	return _get_id($self->get_read_id(0));
#}

sub forward_name {
	my $self = shift;
	return $self->get_read_name(1);
}

sub reverse_name {
	my $self = shift;
	return $self->get_read_name(0);
}

#sub get_read_id{
#	my ( $self, $is_forward ) = @_;
#	my $read_direction = $is_forward ? 1 : 2;
#	#s_$num_$read_direction_sequence.txt
#	my $name =
#	    _get_id($self->{DIR}) . '_' . 's_'
#	  . $self->{NUM} . '_'
#	  . $read_direction
#	  . '_sequence.txt';
#	return $name;
#}

sub get_read_name {
	my ( $self, $is_forward ) = @_;
	my $read_direction = $is_forward ? 1 : 2;
	my $gz = "";
	$gz = '.gz' if $self->{'ARCHIVED'};
	#s_$num_$read_direction_sequence.txt
	my $name =
	    $self->{DIR} . '/' . 's_'
	  . $self->{NUM} . '_'
	  . $read_direction
	  . '_sequence.txt' . $gz;
	return $name;
}

return 1;
