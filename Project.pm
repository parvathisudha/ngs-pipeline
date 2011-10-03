package Project;
use strict;
use Data::Dumper;

#file names for all tasks

sub new {
	my ( $class, $config, $debug ) = @_;
	my $time = get_time($config->{'DIR'});
	if ( $config->{'TIME'} ) {
		$time = $config->{'TIME'};
	}
	my $self = {
		CONFIG => $config->{'PARAMETERS'},
		XML => $config,
		DEBUG  => $debug,
		TIME   => $time,
	};
	bless $self, $class;
	return $self;
}

#folders
sub get_lanes{
	my ($self) = @_;
	my $lanes = $self->{XML}->{DATA}->{LANE};
	return $lanes;
}
sub dir {
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'};
}

sub script_dir {
	my ($self) = @_;
	return $self->dir() . '/tasks';
}

sub output_dir {
	my ($self) = @_;
	return $self->dir() . '/out';
}

sub error_dir {
	my ($self) = @_;
	return $self->dir() . '/errors';
}

sub ids_dir {
	my ($self) = @_;
	return $self->dir() . '/ids';
}

sub make_folder_structure {
	my ($self) = @_;
	$self->mkdir( $self->dir() );
	$self->mkdir( $self->script_dir() );
	$self->mkdir( $self->output_dir() );
	$self->mkdir( $self->error_dir() );
	$self->mkdir( $self->ids_dir() );
	$self->mkdir( $self->tmp_dir() );
}

sub mkdir {
	my ( $self, $directory ) = @_;
	$self->command("mkdir $directory");
}

sub file_prefix {
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->{'CONFIG'}->{'PROJECT'};
}

sub tmp_dir{
	my ($self) = @_;
	return $self->dir() . '/tmp';
}

sub command {
	my ( $self, $command ) = @_;
	if ( $self->{'DEBUG'} ) {
		print "$command\n";
	}
	else {
		system("echo $command");
		system("$command");
	}
}

sub get_time {
	my $dir = shift;
	opendir( DIR, $dir ) or print $!, "\n";
	my @matches = grep( /\d{10}/, readdir(DIR) );
	closedir(DIR);
	if(scalar @matches > 0){
		return $1 if $matches[0] =~ m/(\d{10})/;
	}
	return time;
}

sub read_intervals {
	my ( $self ) = @_;
	my $file = $self->{'CONFIG'}->{'GATKGENOMEBED'};
	my @data;
	open IN, $file or die "Can't open $file\n";
	while (<IN>) {
		chomp;
		next if m/GL/;
		push( @data, "$1" ) if m/(.+)\t(\d+)\t(\d+)/;
	}
	close IN;
	return @data;
}

return 1;