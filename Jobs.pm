use Job;
use Data::Dumper;
#######################################################
{

	package RootJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('root_job');
		$self->virtual(1);
	}
	1;
}

#######################################################
{

	package MergeSamFiles;
	our @ISA = qw( PicardJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('sort_sam');
		my $tmp_dir = $self->project()->dir;
		$self->memory(5);
		my $previous = $self->previous();
		my @input    = $self->input;
		$output = $input[0] . ".joined.bam";
		my $program =
		  $self->picard( "MergeSamFiles.jar", [ @input, "OUTPUT=$output" ] );
		$self->program($program);
		$self->out($output);
	}

	sub input {
		my ( $self, ) = @_;
		return map { "INPUT=" . $_->out } @{ $self->previous };
	}

	sub out {
		my ( $self, $output ) = @_;
		$self->{output_by_type}->{main} = $output if $output;
		return $self->{output_by_type}->{main};
	}
	1;
}

#######################################################
{

	package SaiToBam;
	use Data::Dumper;
	our @ISA = qw( BwaJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('sai_to_bam');
		my $view =
		  $self->project()->{'CONFIG'}->{'SAMTOOLS'} . "/samtools view";
		my $genome   = $self->project()->{'CONFIG'}->{'GENOME'};
		my $previous = $self->previous();
		my $rg       = $self->get_read_group( $self->lane() );
		my $sai      = join( " ", map { $_->out() } @$previous );
		my $fastq = join( " ", map { $_->output_by_type('fastq') } @$previous );
		my $out   = $$previous[0]->out() . ".bam";
		my $command = scalar @$previous == 2 ? 'sampe' : 'samse';
		my $program = Program->new(
			path         => $self->project()->{'CONFIG'}->{'BWA'},
			name         => $command,
			basic_params =>
			  [ $genome, $sai, $fastq, "-r", $rg, "| $view -bS - >", $out ]
		);
		$self->program($program);
		$self->out($out);
	}

	sub get_read_group {
		my ( $self, $lane ) = @_;
		my $rg = "@RG";
		while ( ( $key, $value ) = each %params ) {
			$rg .= "\t$key:$value" if length($key) == 2;
		}
		return $rg;
	}
	1;
}
#######################################################
{

	package GATKJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		my $gatk =
		  $self->project()->{'CONFIG'}->{'GATK'} . "/GenomeAnalysisTK.jar";
		$gatk .= '-R ' . $self->project()->{'CONFIG'}->{'GENOME'};
		$self->{gatk} = $gatk;
		return $self;
	}

	sub gatk {
		my ( $self, %params ) = @_;
		my $gatk = $self->{gatk};
		$gatk = "java -Xmx" . $self->memory() . "g -jar $gatk";
		while ( ( $key, $value ) = each %params ) {
			$gatk .= " $key $value";
		}
		return $gatk;
	}
	1;
}
#######################################################
{

	package PicardJob;
	use Data::Dumper;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		my $program = PicardProgram->new(path=>$self->project()->{'CONFIG'}->{'PICARD'});
		$self->program($program);
		return $self;
	}

	1;
}
#######################################################
{

	package ProcessLane;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	#@Override
	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('process_lane');
		$self->process();
		$self->virtual(1);
	}

	sub process {
		my ( $self, ) = @_;
		my $lane      = $self->{lane};
		my $params    = $self->params();

		#print "PROCESS: ", Dumper $lane;
		my @sai;
		if ( $lane->{FORWARD} ) {
			my $align = Align->new(
				params   => $params,
				previous => [$self],
				lane     => $lane,
				type     => 'FORWARD'
			);
			$align->align_fastq();
			push( @sai, $align );
		}
		if ( $lane->{REVERSE} ) {
			my $align = Align->new(
				params   => $params,
				previous => [$self],
				lane     => $lane,
				type     => 'REVERSE'
			);
			$align->align_fastq();
			push( @sai, $align );
		}
		my $sai_to_bam =
		  SaiToBam->new( lane => $lane, params => $params, previous => \@sai );

		my $sort_sam =
		  SortSam->new( params => $params, previous => [$sai_to_bam] );
		$self->{last_job} = $sort_sam;
	}
	1;
}
#######################################################
{

	package BwaJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('root_job');
		$self->virtual(1);
	}

	sub lane {
		my ( $self, $lane ) = @_;
		$self->{lane} = $lane if $lane;
		return $self->{lane};
	}

	sub genome {
		my ( $self, ) = @_;
		if ( $self->platform eq 'ILLUMINA' ) {
			return $self->{config}->{PARAMETERS}->{GENOME};
		}
		elsif ( $self->platform eq 'Solid' ) {
			return $self->{config}->{PARAMETERS}->{GENOME_COLOR};
		}
		else {
			return undef;
		}
	}

	sub lane {
		my ( $self, $lane ) = @_;
		$self->{lane} = $lane if $lane;
		return $self->{lane};
	}

	sub platform {
		my ( $self, ) = @_;
		return $self->lane->{PL};
	}
	1;
}
#######################################################
{

	package Align;
	our @ISA = qw( BwaJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('align');
		$self->processors( $self->project()->{'CONFIG'}->{'BWA_PROCESSORS'} );
	}

	sub type {
		my ( $self, $type ) = @_;
		$self->{type} = $type if $type;
		return $self->{type};
	}

	sub align_fastq {
		my ( $self, )  = @_;
		my $lane       = $self->lane();
		my $type       = $self->type();
		my $suffix     = $lane->{ID} . ".$type" . '.sai';
		my $prefix     = $self->project()->file_prefix();
		my $in         = $lane->{$type};
		my $out        = "$prefix.$suffix";
		my $proc       = $self->processors();
		my $genome     = $self->genome;
		my $colorspace = "";
		$colorspace = "-c" if $self->platform eq 'Solid';
		my $bwa_priority = 10;
		my $qsub_param   = "-pe mpi $proc -p $bwa_priority";
		my $program      = Program->new(
			prefix       => $self->project()->{'CONFIG'}->{'BWA'},
			name         => 'bwa',
			basic_params => [
				'aln',     '-q 5',  "-t $proc", $colorspace,
				"-f $out", $genome, $in
			]
		);
		$self->program($program);
		$self->qsub_params($qsub_param);
		$self->out($out);
		$self->output_by_type( 'fastq', $in );
	}

	1;
}
#######################################################
{

	package SortSam;
	our @ISA = qw( PicardJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('sort_sam');
		my $tmp_dir = $self->project()->dir;
		$self->memory(5);
		my $previous = $self->previous();
		my $input    = $$previous[0]->out();
		$output = $input . ".sorted.bam";
		$self->program->additional_params(
			[
				"INPUT=$input",      "OUTPUT=$output",
				"CREATE_INDEX=true", "SORT_ORDER=coordinate"
			]
		);
		$self->program->name("SortSam.jar");
		$self->out($output);
	}
	1;
}
#######################################################
