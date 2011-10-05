package Program;

sub new {
	my ( $class, %params ) = @_;
	my $self = {};
	my %all_params = ( %{ $params{params} }, %params );
	while ( my ( $key, $value ) = each %all_params ) {
		$self->{$key} = $value;
	}
	bless $self, $class;
	return $self;
}

sub additional_params {
	my ( $self, $additional_params ) = @_;
	$self->{additional_params} = $additional_params if $additional_params;
	return $self->{additional_params};
}

sub basic_params {
	my ( $self, $basic_params ) = @_;
	$self->{basic_params} = $basic_params if $basic_params;
	return $self->{basic_params};
}

sub name {
	my ( $self, $name ) = @_;
	$self->{name} = $name if $name;
	return $self->{name};
}

sub path {
	my ( $self, $path ) = @_;
	$self->{path} = $path if $path;
	return $self->{path};
}

sub prefix {
	my ( $self, $prefix ) = @_;
	$self->{prefix} = $prefix if $prefix;
	return $self->{prefix};
}

sub memory {
	my ( $self, $memory ) = @_;
	$self->{memory} = $memory if $memory;
	return $self->{memory};
}

sub to_string {
	my ( $self, ) = @_;
	my $path = $self->path;
	$path =~ s/\/$//;
	my $full_program = $path . '/' . $self->name;
	my @all          = (
		$self->prefix, $full_program,
		@{ $self->additional_parameters },
		@{ $self->basic_parameters }
	);
	return join( " ", @all );
}
1;

package JavaProgram;
our @ISA = qw( Program );

sub new {
	my ( $class, %params ) = @_;
	my $self = $class->SUPER::new(%params);
	bless $self, $class;
	return $self;
}

sub prefix {
	my ( $self, ) = @_;
	return "java -jar Xmx" . $self->memory . "g -jar";
}
1;

package PicardProgram;
our @ISA = qw( JavaProgram );

sub new {
	my ( $class, %params ) = @_;
	my $self = $class->SUPER::new(%params);
	bless $self, $class;
	$self->basic_params(
		[
			"TMP_DIR=" . $self->project()->tmp_dir(),
			"VALIDATION_STRINGENCY=SILENT",
			"MAX_RECORDS_IN_RAM=1250000",
		]
	);
	return $self;
}

1;

package GATKProgram;
our @ISA = qw( JavaProgram );

sub new {
	my ( $class, %params ) = @_;
	my $self = $class->SUPER::new(%params);
	bless $self, $class;
	$self->basic_params(
		[
			"TMP_DIR=" . $self->project()->tmp_dir(),
			"VALIDATION_STRINGENCY=SILENT",
			"MAX_RECORDS_IN_RAM=1250000",
		]
	);
	return $self;
}

1;