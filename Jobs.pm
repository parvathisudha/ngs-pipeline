use Job;
use Program;
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
		$self->string_id('RootJob');
		$self->virtual(1);
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
		$self->string_id('SaiToBam');
		my $view =
		  $self->project()->{'CONFIG'}->{'SAMTOOLS'} . "/samtools view";
		my $genome   = $self->project()->{'CONFIG'}->{'GENOME'};
		my $previous = $self->previous();
		my $rg       = $self->get_read_group( $self->lane() );
		my $sai      = join( " ", map { $_->out() } @$previous );
		my $fastq = join( " ", map { $_->output_by_type('fastq') } @$previous );
		my $out   = $$previous[0]->out() . ".bam";
		my $command = scalar @$previous == 2 ? 'sampe' : 'samse';
		$self->program->name('bwa');
		$self->program->basic_params(
			[
				$command, $genome, $sai,              $fastq,
				"-r",     $rg,     "| $view -bS - >", $out
			]
		);
		$self->out($out);
		$self->memory(4);
	}

	sub get_read_group {
		my ( $self, $lane ) = @_;
		my $rg = "@RG";
		while ( ( $key, $value ) = each %$lane ) {
			$rg .= "\\t$key:$value" if length($key) == 2;
		}
		$rg =~ s/^\\t//;
		$rg = '@RG\t' . $rg;
		return "'$rg'";
	}
	1;
}
#######################################################
{

	package GATKJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self =
		  $class->SUPER::new( %params, program => new GATKProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'GATK'} );
		my @b_params = ("-R " . $self->project()->{'CONFIG'}->{'GENOME'},
				"-T " . $self->walker);
		push(@b_params, "-L " . $self->interval) if $self->interval;
		$self->program->basic_params(
			\@b_params
		);
		return $self;
	}

	sub walker {
		my ( $self, $walker ) = @_;
		$self->{walker} = $walker if $walker;
		return $self->{walker};
	}
	sub interval {
		my ( $self, $interval ) = @_;
		$self->{interval} = $interval if $interval;
		return $self->{interval};
	}
	sub string_id {
		my ( $self, ) = @_;
		return $self->walker;
	}
	1;
}
#######################################################
{

	package RealignerTargetCreator;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self =
		  $class->SUPER::new( %params, 'walker' => 'RealignerTargetCreator' );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(1);
		my $previous = $self->previous();
		my $input    = $self->first_previous->out();
		my $output = $input . "." . $self->interval . ".targets";
		my $KGIND = $self->project()->{'CONFIG'}->{'KGIND'};
		$self->program->additional_params( [ "-o $output", "-I $input", "--known $KGIND"] );
		$self->out($output);
	}
	1;
}

#######################################################
{

	package PicardJob;
	use Data::Dumper;
	use Program;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self =
		  $class->SUPER::new( %params, program => new PicardProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'PICARD'} );
		$self->program->tmp_dir( $self->project()->tmp_dir() );
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
		$self->string_id('ProcessLane');
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
		$self->program->path( $self->project()->{'CONFIG'}->{'BWA'} );
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
		$self->string_id('Align');
		$self->processors( $self->project()->{'CONFIG'}->{'BWA_PROCESSORS'} );
		$self->memory(5);
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
		$self->program->name(bwa);
		my $illumina_fastq_qualities_flag = "";
		$illumina_fastq_qualities_flag = "-I" if $self->lane->{ILLUMINA_FASTQ};
		$self->program->basic_params(
			[ 'aln', '-q 5', "-t $proc", $illumina_fastq_qualities_flag ,$colorspace, "-f $out", $genome, $in ]
		);
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
		$self->string_id('SortSam');
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
{

	package MergeSamFiles;
	our @ISA = qw( PicardJob );
	use Data::Dumper;

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('MergeSamFiles');
		$self->memory(5);
		my $previous = $self->previous();
		my @input    = map { "INPUT=" . $_->out } @$previous;
		my $input    = join( " ", @input );
		my $output = $self->out;
		$self->program->additional_params(
			[
				"$input",            "OUTPUT=$output",
				"CREATE_INDEX=true", "SORT_ORDER=coordinate"
			]
		);
		$self->program->name("MergeSamFiles.jar");
	}
	1;
}
#######################################################
{

	package MarkDuplicates;
	our @ISA = qw( PicardJob );
	use Data::Dumper;

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('MarkDuplicates');
		$self->memory(5);
		my $previous = $self->previous();
		my $input = $self->first_previous->out;
		my $output   = $input . ".dedup.bam";
		my $metrics  = $input . ".metrics.txt";
		$self->program->additional_params(
			[
				"INPUT=$input" , "OUTPUT=$output",
				"CREATE_INDEX=true", "SORT_ORDER=coordinate",
				"METRICS_FILE=$metrics",
			]
		);
		$self->program->name("MarkDuplicates.jar");
		$self->out($output);
		$self->output_by_type( "metrics", $metrics );
	}
	1;
}
#######################################################
