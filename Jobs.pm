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

	package SaiToBam;
	use Data::Dumper;
	our @ISA = qw( Job );
	
	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->string_id('sai_to_bam');
		my $bwa_path = $self->project()->{'CONFIG'}->{'BWA'};
		my $view     =
		  $self->project()->{'CONFIG'}->{'SAMTOOLS'} . "/samtools view";
		my $genome  = $self->project()->{'CONFIG'}->{'GENOME'};
		my $previous = $self->previous();
		my $rg      = $self->get_read_group( $$previous[0]->lane() );
		my $sai     = join( " ", map { $_->out() } @$previous );
		my $fastq   = join( " ", map { $_->output_by_type('fastq') } @$previous );
		my $out     = $$previous[0]->out() . ".bam";
		my $command = scalar @$previous == 2 ? 'sampe' : 'samse';
		$program =
		  "$bwa_path/$command $genome $sai $fastq -r $rg | $view -bS - > $out";
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
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		my $picard = $self->project()->{'CONFIG'}->{'PICARD'} . "/";
		$self->{picard} = $picard;
		return $self;
	}

	sub picard {
		my ( $self, %params ) = @_;
		my $picard = $self->{picard};
		$picard = "java -Xmx" . $self->memory() . "g -jar $picard";
		$picard .= $params{jar};
		delete $params{jar};
		my %basic_params = {
			TMP_DIR               => $self->project()->tmp_dir(),
			VALIDATION_STRINGENCY => SILENT,
			MAX_RECORDS_IN_RAM    => 1250000,
		};
		my %all_params = ( %basic_params, %params );
		while ( ( $key, $value ) = each %all_params ) {
			$picard .= " $key=$value";
		}
		return $picard;
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
		my $job_factory = $self->job_factory();
		if ( $lane->{PL} eq 'ILLUMINA' ) {
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
			  SaiToBam->new( params => $params, previous =>  \@sai );
			
			my $sort_sam =
			  SortSam->new( params => $params, previous => [$sai_to_bam] );
			$self->{last_job} = $sort_sam;
		}
	}

	#@Override
	sub output_files {
		my ( $self, ) = @_;
	}

	#@Override
	sub program {
		my ( $self, ) = @_;
	}
	1;
}
#######################################################
{

	package Align;
	our @ISA = qw( Job );

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

	sub lane {
		my ( $self, $lane ) = @_;
		$self->{lane} = $lane if $lane;
		return $self->{lane};
	}
	sub type {
		my ( $self, $type ) = @_;
		$self->{type} = $type if $type;
		return $self->{type};
	}
	sub align_fastq {
		my ( $self,) = @_;
		my $lane = $self->lane();
		my $type = $self->type();
		my $suffix       = $lane->{ID} . ".$type" . '.sai';
		my $prefix       = $self->project()->file_prefix();
		my $in           = $lane->{$type};
		my $out          = "$prefix.$suffix";
		my $bwa          = $self->project()->{'CONFIG'}->{'BWA'} . "/bwa";
		my $proc         = $self->processors();
		my $align        = "$bwa aln -t $proc";
		my $program      = "$align -f $out $genome $in";
		my $bwa_priority = 10;
		my $qsub_param   = "-pe mpi $proc -p $bwa_priority";
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
		my $input = $$previous[0]->out();
		$output = $input . "sorted.bam";
		my $program = $self->picard(
			jar          => SortSam . jar,
			INPUT        => $input,
			OUTPUT       => $output,
			CREATE_INDEX => true,
			SORT_ORDER   => coordinate,
		);
		$self->program($program);
		$self->out($output);
	}
	1;
}
#######################################################
