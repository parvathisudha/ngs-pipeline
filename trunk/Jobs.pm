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
		my $view =
		  $self->project()->{'CONFIG'}->{'SAMTOOLS'} . "/samtools view";
		my $genome   = $self->first_previous->output_by_type('genome');
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

	package SnpEff;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, program => new JavaProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'SNPEFF'} );
		$self->program->name("snpEff.jar");
		$self->memory(8);
		my $input = $self->first_previous->output_by_type('vcf');
		my $html = $input . ".html";
		my $snpeff_genome = $self->project->{'CONFIG'}->{SNPEFF_GENOME};
		my $output = $input . ".eff.vcf";
		$self->program->additional_params(
			[ "-onlyCoding", "-stats $html", "-o vcf", $snpeff_genome, "$input > $output" ] );
		$self->out($output);
		$self->output_by_type( 'vcf', $output );
		return $self;
	}

	1;
}
#######################################################
{

	package GATKJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, program => new GATKProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'GATK'} );
		$self->walker($class);
		my @b_params = (
			"-R " . $self->project()->{'CONFIG'}->{'GENOME'},
			"-T " . $self->walker
		);
		push( @b_params, "-L " . $self->interval ) if $self->interval;
		$self->program->basic_params( \@b_params );
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
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $previous = $self->previous();
		my $input    = $self->first_previous->out();
		my $output   = $input . "." . $self->interval . ".targets";
		my $KGIND    = $self->project()->{'CONFIG'}->{'KGIND'};
		$self->program->additional_params(
			[ "-o $output", "-I $input", "--known $KGIND" ] );
		$self->out($output);
		$self->output_by_type( 'bam', $input );
	}
	1;
}
#######################################################
{

	package IndelRealigner;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $previous = $self->previous();
		my $input    = $self->first_previous->output_by_type('bam');
		my $targets  = $self->first_previous->out;
		my $output   = $input . "." . $self->interval . ".realigned.bam";
		my $KGIND    = $self->project()->{'CONFIG'}->{'KGIND'};
		$self->program->additional_params(
			[
				"-o $output",
				"-I $input",
				"--known $KGIND",
				"-targetIntervals $targets"
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $input );
	}
	1;
}
#######################################################
{

	package CountCovariates;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input  = $self->first_previous->output_by_type('bam');
		my $output = $input . ".covarities";
		my $dbSNP  = $self->project()->{'CONFIG'}->{'DBSNP'};
		$self->program->additional_params(
			[
				"-I $input",
				"-cov ReadGroupCovariate",
				"-cov QualityScoreCovariate",
				"-cov CycleCovariate",
				"-cov DinucCovariate",
				"-recalFile $output",
				"-knownSites $dbSNP",
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $input );
	}
	1;
}
#######################################################
{

	package TableRecalibration;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input      = $self->first_previous->output_by_type('bam');
		my $recal_file = $self->first_previous->out;
		my $output     = $input . ".recalibrated.bam";
		my $dbSNP      = $self->project()->{'CONFIG'}->{'DBSNP'};
		$self->program->additional_params(
			[
				"-l INFO",
				"-o $output",
				"-I $input",
				"-recalFile $recal_file",
				"-knownSites $dbSNP",
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $output );
	}
	1;
}
#######################################################
{

	package SelectVariants;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}
	sub sample{
		my ($self) = @_;
		return $self->{sample};
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input      = $self->in ? $self->in : $self->first_previous->output_by_type('vcf');
		$self->program->additional_params(
			[
				"-o " . $self->out,
				"--variant $input",
			]
		);
		$self->output_by_type( 'vcf', $self->out );
	}
	1;
}
#######################################################
{

	package UnifiedGenotyper;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub variation_type {
		my ( $self, ) = @_;
		return $self->{variation_type};
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input           = $self->first_previous->output_by_type('bam');
		my $recal_file      = $self->first_previous->out;
		my $output          = $input . "." . $self->variation_type . ".vcf";
		my $dbSNP           = $self->project()->{'CONFIG'}->{'DBSNP'};
		my $stand_call_conf =
		  $self->project()->{'CONFIG'}->{'GATK_stand_call_conf'};
		my $stand_emit_conf =
		  $self->project()->{'CONFIG'}->{'GATK_stand_emit_conf'};
		my @annotations = qw/
		  ChromosomeCounts
		  IndelType
		  HardyWeinberg
		  SpanningDeletions
		  GLstats
		  NBaseCount
		  AlleleBalance
		  MappingQualityZero
		  BaseCounts
		  LowMQ
		  RMSMappingQuality
		  HaplotypeScore
		  TechnologyComposition
		  SampleList
		  FisherStrand
		  DepthOfCoverage
		  HomopolymerRun
		  MappingQualityZeroFraction
		  GCContent
		  MappingQualityRankSumTest
		  ReadPosRankSumTest
		  BaseQualityRankSumTest
		  QualByDepth
		  SBByDepth
		  ReadDepthAndAllelicFractionBySample
		  AlleleBalanceBySample
		  DepthPerAlleleBySample
		  MappingQualityZeroBySample
		  /;
		my @formatted_annotations = map { "--annotation $_" } @annotations;
		my $genotype_likelihoods_model = $self->variation_type;
		$self->program->additional_params(
			[
				"-o $output",
				"-I $input",
				"--dbsnp $dbSNP",
				"-stand_call_conf $stand_call_conf",
				"-stand_emit_conf $stand_emit_conf",
				"-dcov 80 -U",
				"--genotype_likelihoods_model $genotype_likelihoods_model",
				@formatted_annotations,
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $input );
		$self->output_by_type( 'vcf', $output );
	}
	1;
}
#######################################################
{

	package VariantAnnotator;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input      = $self->first_previous->output_by_type('bam');
		my $variations = $self->first_previous->output_by_type('vcf');
		my $output     = $variations . ".annotated.vcf";
		my $bam = "";
		$bam = "-I $input" if $input;
		$self->program->additional_params(
			[
				"-o $output",
				$bam,
				"--variant $variations",
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $input );
		$self->output_by_type( 'vcf', $output );
	}
	1;
}
#######################################################
{

	package CombineVariants;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $sample   = $self->project->{CONFIG}->{SAMPLE_NAME};
		my @previous = @{ $self->previous };
		my @vcfs     =
		  map { "--variant:$sample " . $_->output_by_type('vcf') } @previous;
		my $input = join( " ", @vcfs );
		my $output = $self->out;
		$self->program->additional_params(
			[ $input, "-o $output", "--assumeIdenticalSamples", ] );
		$self->out($output);
		$self->output_by_type( 'vcf', $output );
	}
	1;
}
#######################################################
{

	package VariantRecalibrator;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(4);
		my $input         = $self->first_previous->output_by_type('vcf');
		my $recal_file    = $input . ".recal";
		my $tranches_file = $input . ".tranches";
		my $rscript_file  = $input . ".plots.R";
		my $resources     = $self->project->{CONFIG}->{GATK} . "/resources/";
		$self->program->additional_params(
			[
				"-input $input",
				"--recal_file $recal_file",
				"--tranches_file $tranches_file",
				"--rscript_file $rscript_file",
				"--path_to_resources $resources",
			]
		);
		$self->out($tranches_file);
		$self->output_by_type( 'vcf',           $input );
		$self->output_by_type( 'recal_file',    $recal_file );
		$self->output_by_type( 'tranches_file', $tranches_file );
		$self->output_by_type( 'rscript_file',  $rscript_file );
	}
	1;
}
#######################################################
{

	package ApplyRecalibration;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(6);
		my $input  = $self->first_previous->output_by_type('vcf');
		my $output = $input . ".recalibrated.vcf";
		$self->program->additional_params(
			[
				"-input $input",
				"--ts_filter_level 99.0",
				"--recal_file "
				  . $self->first_previous->output_by_type('recal_file'),
				"--tranches_file "
				  . $self->first_previous->output_by_type('tranches_file'),
				"-o $output"
			]
		);
		$self->out($output);
		$self->output_by_type( 'vcf',    $output );
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
		$self->program->name( $class . ".jar" );
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
			[
				'aln', '-q 5', "-t $proc", $illumina_fastq_qualities_flag,
				$colorspace, "-f $out", $genome, $in
			]
		);
		$self->qsub_params($qsub_param);
		$self->out($out);
		$self->output_by_type( 'fastq',  $in );
		$self->output_by_type( 'genome', $genome );
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
		$self->out($output);
		$self->output_by_type( 'bam', $output );
	}
	1;
}
#######################################################
{

	package BuildBamIndex;
	our @ISA = qw( PicardJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		my $tmp_dir = $self->project()->dir;
		$self->memory(5);
		my $input  = $self->first_previous->out();
		my $output = $input;
		$output =~ s/bam$/bai/;
		$self->program->additional_params(
			[ "INPUT=$input", "OUTPUT=$output", ] );
		$self->out($output);
		$self->output_by_type( 'bam', $input );
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
		$self->memory(5);
		my $previous = $self->previous();
		my @input    = map { "INPUT=" . $_->out } @$previous;
		my $input    = join( " ", @input );
		my $output   = $self->out;
		$self->program->additional_params(
			[
				"$input",            "OUTPUT=$output",
				"CREATE_INDEX=true", "SORT_ORDER=coordinate"
			]
		);
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
		$self->memory(5);
		my $previous = $self->previous();
		my $input    = $self->first_previous->out;
		my $output   = $input . ".dedup.bam";
		my $metrics  = $input . ".metrics.txt";
		$self->program->additional_params(
			[
				"INPUT=$input",      "OUTPUT=$output",
				"CREATE_INDEX=true", "METRICS_FILE=$metrics",
			]
		);
		$self->out($output);
		$self->output_by_type( "metrics", $metrics );
	}
	1;
}
#######################################################
