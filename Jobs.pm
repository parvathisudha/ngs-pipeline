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
	our @ISA = qw( Aligner );

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
		my $rg = '@RG';
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

	package IGVTools;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, program => new JavaProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'IGVTOOLS'} );
		$self->program->name("igvtools.jar");
		$self->memory(2);
		return $self;
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
		my $config =
		  $self->project()->{'CONFIG'}->{'SNPEFF'} . "/snpEff.config";
		$self->program->name("snpEff.jar");
		$self->memory(8);
		my $input         = $self->first_previous->output_by_type('vcf');
		my $html          = $input . ".html";
		my $snpeff_genome = $self->project->{'CONFIG'}->{SNPEFF_GENOME};
		my $output        = $input . ".eff.vcf";
		$self->program->additional_params(
			[
				"-config $config",
				"-onlyCoding true",
				"-stats $html",
				"-o vcf",
				$snpeff_genome,
				"$input -reg base > $output"
			]
		);
		$self->out($output);
		$self->output_by_type( 'vcf', $output );
		return $self;
	}

	1;
}
#######################################################
{

	package FilterLowQual;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("grep");
		$self->program->path("/bin");
		$self->memory(1);
		my $input  = $self->first_previous->output_by_type('vcf');
		my $output = $input . ".filtered.vcf";
		$self->program->additional_params( ["-v LowQual $input > $output"] );
		$self->out($output);
		$self->output_by_type( 'vcf', $output );
		return $self;
	}
	1;
}
#######################################################
{

	package Grep;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("grep");
		$self->program->path("/bin");
		$self->memory(1);
		return $self;
	}
	1;
}
#######################################################
{

	package Gunzip;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("gunzip");
		$self->program->path("/bin");
		$self->memory(1);
		my $input  = $self->in;
		my $output = $self->out;
		$self->program->additional_params( ["-cd $input > $output"] );
		return $self;
	}
	1;
}
#######################################################
{

	package Substitute;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("sed");
		$self->program->path("/bin");
		die "The input file is not specified!\n",  unless $self->in;
		die "The output file is not specified!\n", unless $self->out;
		$self->memory(1);
		$self->program->additional_params(
			[
				"\'s/" . $self->{from} . '/' . $self->{to} . '/\'',
				$self->in, ">", $self->out
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package Ln;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("ln");
		$self->program->path("/bin");
		die "The input file is not specified!\n",  unless $self->in;
		die "The output file is not specified!\n", unless $self->out;
		$self->memory(1);
		$self->output_by_type( 'link', $self->out );
		$self->program->additional_params( [ "-s", $self->in, $self->out ] );
		return $self;
	}
	1;
}
#######################################################
{

	package Pindel2Vcf;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("pindel2vcf");
		$self->program->path( $self->project()->{'CONFIG'}->{'PINDEL'} );
		$self->memory(1);
		$self->output_by_type( 'vcf', $self->out );
		$self->program->additional_params(
			[
				"-G",    #only for the private version of pindel2vcf
				"--reference",
				$self->project()->{'CONFIG'}->{'GENOME'},
				"--reference_name",
				$self->project()->{'CONFIG'}->{'GENOME_NAME'},
				"--reference_date",
				$self->project()->{'CONFIG'}->{'DATE'},
				"--pindel_output",
				$self->in,
				"--vcf",
				$self->out,
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package Pindel;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("pindel");
		$self->program->path( $self->project()->{'CONFIG'}->{'PINDEL'} );
		$self->memory(8);
		my $prefix     = $self->project->file_prefix() . ".pindel";
		my $out        = $prefix;
		my $threads    = 4;
		my $qsub_param = "-pe mpi $threads";
		$self->qsub_params($qsub_param);
		$self->out($out);
		$self->output_by_type( 'SV_D',   $prefix . "_D" );
		$self->output_by_type( 'SV_INV', $prefix . "_INV" );
		$self->output_by_type( 'SV_LI',  $prefix . "_LI" );
		$self->output_by_type( 'SV_SI',  $prefix . "_SI" );
		$self->output_by_type( 'SV_TD',  $prefix . "_TD" );
		$self->output_by_type( 'SV_BP',  $prefix . "_BP" );
		$self->program->additional_params(
			[
				"--fasta",
				$self->project()->{'CONFIG'}->{'GENOME'},
				"--config-file",
				$self->first_previous->output_by_type('cfg'),
				"--chromosome ALL",
				"--output-prefix $prefix",
				"--number_of_threads $threads",
				"> $out"
			]
		);
		return $self;
	}

	sub variation_files {
		my ($self) = @_;
		my @result;
		while ( ( $type, $file ) = each %{ $self->{output_by_type} } ) {
			push( @result, $file ) if $type =~ m/SV/;
		}
		return \@result;
	}

	sub VEP_compatible_files {
		my ($self) = @_;
		my @result;
		while ( ( $type, $file ) = each %{ $self->{output_by_type} } ) {
			push( @result, $file )
			  if ( $type =~ m/SV_D/
				|| $type =~ m/SV_SI/
				|| $type =~ m/SV_INV/ );    #
		}
		return \@result;
	}

	sub SnpEff_compatible_files {
		my ($self) = @_;
		my @result;
		while ( ( $type, $file ) = each %{ $self->{output_by_type} } ) {
			push( @result, $file )
			  if ( $type =~ m/SV_D/ || $type =~ m/SV_LI/ || $type =~ m/SV_SI/ )
			  ;                             #
		}
		return \@result;
	}
	1;
}
#######################################################
{

	package Cat;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("cat");
		$self->program->path("/bin");
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		return $self;
	}
	1;
}
#######################################################
{

	package Cut;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("cut");
		$self->program->path("/bin");
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		return $self;
	}
	1;
}
#######################################################
{

	package Bgzip;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("bgzip");
		$self->program->path( $self->project()->{'CONFIG'}->{'TABIX'} );
		$self->memory(1);
		my $vcf = $self->first_previous->output_by_type('vcf');
		my $gz  = "$vcf.gz";
		$self->output_by_type( 'gz', $gz );
		$self->out($gz);
		$self->output_by_type( 'vcf', $vcf );
		$self->program->additional_params( ["-c $vcf > $gz"] );
		return $self;
	}
	1;
}
#######################################################
{

	package GenomeCoverageBed;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("genomeCoverageBed");
		$self->program->path( $self->project()->{'CONFIG'}->{'BEDTOOLS'} );
		$self->memory(1);
		my $bam = $self->first_previous->output_by_type('bam');
		my $cov = "$bam.cov";
		$self->output_by_type( 'cov', $cov );
		$self->output_by_type( 'bam', $bam );
		$self->out($cov);
		my $genome = $self->project()->{'CONFIG'}->{'BEDGENOME'};
		$self->program->additional_params( ["ibam $bam -g $genome"] );
		return $self;
	}
	1;
}
#######################################################
{

	package Tabix;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("tabix");
		$self->program->path( $self->project()->{'CONFIG'}->{'TABIX'} );
		$self->memory(1);
		my $gz  = $self->first_previous->output_by_type('gz');
		my $tbi = "$gz.tbi";
		$self->output_by_type( 'gz',  $gz );
		$self->output_by_type( 'tbi', $tbi );
		$self->out($tbi);
		$self->program->additional_params( ["-p vcf $gz"] );
		return $self;
	}
	1;
}
#######################################################
{

	package FREEC;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("freec");
		$self->program->path(
			$self->project()->{'CONFIG'}->{'FREEC'}->{'PATH'} );
		$self->memory(1);
		my $conf        = $self->first_previous->output_by_type('txt');
		my $name_prefix =
		  $self->first_previous->first_previous->output_by_type('bam');

		#		my $cpn = $name_prefix . "_sample.cpn";
		my $cnv = $name_prefix . "_CNVs";
		$self->output_by_type( 'cnv',   $cnv );
		$self->output_by_type( 'ratio', $name_prefix . "_ratio.txt" );
		$self->out($cnv);
		$self->program->additional_params( ["-conf $conf"] );
		return $self;
	}
	1;
}
#######################################################
{

	package FREEC2Circos;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("freec2circos.pl");
		$self->program->path(
			$self->project()->{'CONFIG'}->{'FREEC'}->{'SCRIPTS'} );
		$self->memory(1);
		my $out =
		  $self->first_previous()->output_by_type('ratio') . ".circos.hist";
		$self->output_by_type( 'conf', $out );
		$self->out($out);
		$self->program->additional_params(
			[
				"-f", $self->first_previous()->output_by_type('ratio'),
				"-p", $self->project()->{'CONFIG'}->{'FREEC'}->{'ploidy'},
				"> $out"
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package VcfSorter;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("vcfsorter.pl");
		$self->program->path( $self->project()->{'CONFIG'}->{'VCF_SORTER'} );
		$self->memory(8);
		my $vcf    = $self->first_previous->output_by_type('vcf');
		my $dict   = $self->project()->{'CONFIG'}->{'GENOME_DICT'};
		my $sorted = "$vcf.s.vcf";
		$self->output_by_type( 'vcf', $sorted );
		$self->output_by_type( 'idx', $sorted . '.idx' );
		$self->out($sorted);
		$self->program->additional_params( ["$dict $vcf > $sorted"] );
		return $self;
	}
	1;
}

#######################################################
{

	package VcfIntersectionToReport;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("vcf_intersection_to_report.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->program->additional_params( [ "--out", $self->out ] );
		$self->memory(1);
		return $self;
	}
	1;
}

#######################################################
{

	package FilterFreq;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("filter_vcf.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $vcf  = $self->first_previous->output_by_type('vcf');
		my $rare = "$vcf.rare.vcf";
		$self->output_by_type( 'vcf', $rare );
		$self->out($rare);
		$self->program->additional_params( ["$vcf $rare"] );
		return $self;
	}
	1;
}
#######################################################
{

	package IntersectVcfBed;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("intersect.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $vcf  = $self->first_previous->output_by_type('vcf');
		my $rare = $self->out;
		$self->output_by_type( 'vcf', $rare );

		$self->program->additional_params(
			[
				"--in $vcf", "--bedtools",
				$self->project()->{'CONFIG'}->{'BEDTOOLS'},
				"--bed", $self->{bed}, "--out", $rare
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package PerlJob;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, program => new PerlProgram() );
		bless $self, $class;

		#my $perl_lib = join (" ", map {"-I $_"} @{$self->{perl_lib}});
		my $perl_lib = $self->project()->{'CONFIG'}->{'PERLLIB'};
		$self->program->prefix("perl -I $perl_lib");
		return $self;
	}
	1;
}
#######################################################
{

	package PindelConfig;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("pindel_config.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $in  = $self->first_previous->output_by_type('txt');
		my $out = "$in.pindel";
		$self->output_by_type( 'cfg', $out );
		$self->out($out);
		$self->program->additional_params( ["--in $in --out $out"] );
		return $self;
	}
	1;
}
#######################################################
{

	package GrepTxt;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("grep_txt.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $txt = $self->first_previous->output_by_type('txt');
		my $grep = $self->out ? $self->out : "$txt.grep.txt";
		$self->output_by_type( 'txt', $grep );
		$self->out($grep);
		$self->program->additional_params( ["--in $txt --out $grep"] );
		return $self;
	}
	1;
}
#######################################################
{

	package GrepVcf;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("grep_vcf.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $vcf = $self->first_previous->output_by_type('vcf');
		my $grep = $self->out ? $self->out : "$vcf.grep.vcf";
		$self->output_by_type( 'vcf', $grep );
		$self->out($grep);
		$self->program->additional_params( ["--in $vcf --out $grep"] );
		return $self;
	}
	1;
}
#######################################################
{

	package VEP;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("variant_effect_predictor.pl");
		$self->program->path( $self->project()->{'CONFIG'}->{'VEP'} );
		$self->memory(4);

		my $vcf = $self->first_previous->output_by_type('vcf');
		my $vep = $self->out;
		$self->output_by_type( 'vcf', $vep );
		my $proc = 6;

		#		my $qsub_param   = "-pe mpi $proc";
		#		$self->qsub_params($qsub_param);
		$self->program->additional_params(
			[
				"--input_file $vcf",
				"--format vcf",
				"--cache",
				"--vcf",
				"--buffer_size 30000",
				"--force_overwrite",
				"--dir " . $self->project()->{'CONFIG'}->{'VEPCACHE'},
				"--output_file $vep",

				#				"--fork $proc",
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package CodingReport;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("coding_report.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $vcf   = $self->first_previous->output_by_type('vcf');
		my $table = $self->out;
		$self->output_by_type( 'txt', $table );
		$self->program->additional_params(
			[
				"--in $vcf",
				"--out $table",
				"--annotation_types_file",
				$self->project()->install_dir
				  . "/accessory/annotation_types.txt"
			]
		);
		return $self;
	}
	1;
}
#######################################################
{

	package VcfToReport;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("vcf_to_report.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $vcf   = $self->first_previous->output_by_type('vcf');
		my $table = $self->out;
		$self->output_by_type( 'xls', $table );
		$self->program->additional_params( [ "--in $vcf", "--out $table", ] );
		return $self;
	}
	1;
}
#######################################################
{

	package WriteFreecConfig;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("write_freec_config.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		return $self;
	}
	1;
}
#######################################################
{

	package AddLoci;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->program->name("add_loci.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $in     = $self->first_previous->output_by_type('txt');
		my $result = $in . '.loci.txt';
		$self->out($result);
		$self->output_by_type( 'txt', $self->out );
		$self->program->additional_params( [ "--in $in", "--out $result", ] );
		return $self;
	}
	1;
}
#######################################################
{

	package ReformatRegulation;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("reformat_regulation_info.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		$self->program->basic_params( [ ">", $self->out, ] );
		return $self;
	}
	1;
}
#######################################################
{

	package AnnotateProteins;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("gene_to_annotation.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		$self->program->basic_params( [ ">", $self->out, ] );
		return $self;
	}
	1;
}
#######################################################
{

	package JoinTabular;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("join_tabular_files.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		$self->output_by_type( 'txt', $self->out );
		$self->program->basic_params( [ ">", $self->out, ] );
		return $self;
	}
	1;
}
#######################################################
{

	package NormalizeTelomeres;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("normalize_telomers.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $out = $self->first_previous->out . ".norm";
		$self->out($out);
		$self->program->additional_params(
			[
				"--normalization",
				$self->project()->install_dir
				  . "/accessory/normalization.tel.distribution",
				"--input",
				$self->first_previous->out,
				"--result $out",
			]
		);
		return $self;
	}
	1;
}

#######################################################
{

	package CalcTelomeres;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("calc_telomeres.pl");
		$self->program->path( $self->project()->install_dir . "/accessory" );
		$self->memory(1);
		my $lane         = $self->{lane};
		my $id           = $lane->{ID};
		my $file_name    = $self->project()->file_prefix() . "." . $id . ".tel";
		my $distribution = $file_name . ".result";
		my $samtools     = $self->project()->{CONFIG}->{SAMTOOLS} . "/samtools";
		my $fastq_quality_filter =
		  $self->project()->{CONFIG}->{FASTQ_QUALITY_FILTER};
		$self->output_by_type( 'distribution', $distribution );
		$self->out($distribution);
		$self->program->additional_params(
			[
				"--reads_limit 100000000",
				"--distribution $distribution",
				"--samtools $samtools",
				"--fastq_quality_filter $fastq_quality_filter",
			]
		);

		if ( $lane->{FORWARD} ) {
			$self->program->additional_params( [ "--file", $lane->{FORWARD} ] );
		}
		if ( $lane->{REVERSE} ) {
			$self->program->additional_params( [ "--file", $lane->{REVERSE} ] );
		}
		if ( $lane->{BAM1} ) {
			$self->program->additional_params( [ "--file", $lane->{BAM1} ] );
		}

		return $self;
	}
	1;
}

#######################################################
{

	package CNVSeq;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("cnv-seq.pl");
		$self->program->path( $self->project()->{'CONFIG'}->{'CNVSEQ'} );
		$self->memory(4);
		$self->output_by_type( 'txt', $self->out );
		return $self;
	}
	1;
}

#######################################################
{

	package Bam2cfg;
	our @ISA = qw( PerlJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("bam2cfg.pl");
		$self->program->path( $self->project()->{'CONFIG'}->{'BREAKDANCER'} );
		$self->memory(1);
		my $input = $self->first_previous->output_by_type('bam');
		$self->output_by_type( 'bam', $input );
		my $output = $input . '.cfg';
		$self->output_by_type( 'cfg', $output );
		$self->out($output);
		$self->program->additional_params( ["$input > $output"] );
		return $self;
	}
	1;
}
#######################################################
{

	package BreakdancerMax;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->program->name("breakdancer-max");
		$self->program->path( $self->project()->{'CONFIG'}->{'BREAKDANCER'} );
		$self->memory(5);
		my $bam = $self->first_previous->output_by_type('bam');
		my $cfg = $self->first_previous->output_by_type('cfg');

		my $max   = $cfg . '.max';
		my $fastq = $cfg . '.fastq';
		my $bed   = $cfg . '.bed';

		$self->output_by_type( 'bam',   $bam );
		$self->output_by_type( 'cfg',   $cfg );
		$self->output_by_type( 'max',   $max );
		$self->output_by_type( 'fastq', $fastq );
		$self->output_by_type( 'bed',   $bed );
		$self->out($max);

		$self->program->additional_params( ["-d $fastq -g $bed $cfg > $max"] );

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

	sub idx_from_vcf {
		my ( $self, $vcf ) = @_;
		return "$vcf.idx";
	}

	1;
}
#######################################################
{

	package DepthOfCoverage;
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
		my $input         = $self->first_previous->output_by_type('bam');
		my $output_prefix = $input . ".cov";
		$self->program->additional_params(
			[ "-o $output_prefix", "-I $input", ] );
		$self->out( $output_prefix . '.sample_summary' );
		$self->output_by_type( 'sample_summary', $self->out );
		$self->output_by_type( 'sample_statistics',
			$output_prefix . '.sample_statistics' );
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
		my $output   = $input . "." . $self->interval . ".intervals";
		my $KGIND    = $self->project()->{'CONFIG'}->{'KGIND'};
		$self->program->additional_params(
			[ "-o $output", "-I $input", "-known $KGIND" ] );
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
				"-known $KGIND",
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
				"--solid_nocall_strategy PURGE_READ",
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
				"-l INFO", "-o $output", "-I $input",
				"-recalFile $recal_file",
				"--solid_nocall_strategy PURGE_READ",
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $output );
	}
	1;
}
#######################################################
{

	package VariantsToTable;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		$self->memory(2);
		my $input =
		  $self->in ? $self->in : $self->first_previous->output_by_type('vcf');
		my $out = $self->out ? $self->out : "$input.txt";
		$self->program->basic_params(
			[ "-o " . $out, "-V $input", "--allowMissingData", ] );
		$self->output_by_type( 'txt', $out );
		$self->out($out);
		return $self;
	}

	sub sample {
		my ($self) = @_;
		return $self->{sample};
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

	sub sample {
		my ($self) = @_;
		return $self->{sample};
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(2);
		my $input =
		  $self->in ? $self->in : $self->first_previous->output_by_type('vcf');
		$self->program->additional_params(
			[ "-o " . $self->out, "--variant $input", ] );
		$self->output_by_type( 'vcf', $self->out );
	}
	1;
}
#######################################################
{

	package VariantFiltration;
	our @ISA = qw( SelectVariants );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	1;
}
#######################################################
{

	package ReadBackedPhasing;
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
		my $bam =
		    $self->bam
		  ? $self->bam
		  : $self->first_previous->output_by_type('bam');
		my $input  = $self->first_previous->output_by_type('vcf');
		my $output = $input . ".phased.vcf";
		$self->program->additional_params(
			[
				"-I $bam",
				"--variant $input",

				#				"-L $input",
				"-o $output",
				"-enableMergeToMNP",
				"--maxGenomicDistanceForMNP 2",
			]
		);
		$self->out($output);
		$self->output_by_type( 'vcf', $output );
		$self->output_by_type( 'bam', $bam );
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
		$self->memory(8);
		my $input           = $self->first_previous->output_by_type('bam');
		my $recal_file      = $self->first_previous->out;
		my $output          = $input . "." . $self->variation_type . ".vcf";
		my $dbSNP           = $self->project()->{'CONFIG'}->{'DBSNP'};
		my $stand_call_conf =
		  $self->project()->{'CONFIG'}->{'GATK_stand_call_conf'};
		my $stand_emit_conf =
		  $self->project()->{'CONFIG'}->{'GATK_stand_emit_conf'};

		#		my @annotations = qw/
		#		  ChromosomeCounts
		#		  IndelType
		#		  HardyWeinberg
		#		  SpanningDeletions
		#		  GLstats
		#		  NBaseCount
		#		  AlleleBalance
		#		  MappingQualityZero
		#		  BaseCounts
		#		  LowMQ
		#		  RMSMappingQuality
		#		  HaplotypeScore
		#		  TechnologyComposition
		#		  SampleList
		#		  FisherStrand
		#		  DepthOfCoverage
		#		  HomopolymerRun
		#		  MappingQualityZeroFraction
		#		  GCContent
		#		  MappingQualityRankSumTest
		#		  ReadPosRankSumTest
		#		  BaseQualityRankSumTest
		#		  QualByDepth
		#		  SBByDepth
		#		  ReadDepthAndAllelicFractionBySample
		#		  AlleleBalanceBySample
		#		  DepthPerAlleleBySample
		#		  MappingQualityZeroBySample
		#		  /;
		#		my @formatted_annotations = map { "--annotation $_" } @annotations;
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

				#				@formatted_annotations,
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $input );
		$self->output_by_type( 'vcf', $output );
		$self->output_by_type( 'idx', $self->idx_from_vcf($output) );

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
		$self->memory(30);
		my $input      = $self->first_previous->output_by_type('bam');
		my $variations = $self->first_previous->output_by_type('vcf');
		my $output = $self->out ? $self->out : $variations . ".annotated.vcf";

		my $bam = "";
		$bam = "-I $input" if $input;
		$self->program->additional_params(
			[
				"-o $output",            $bam,
				"--variant $variations", "--requireStrictAlleleMatch"
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

	package LeftAlignVariants;
	our @ISA = qw( GATKJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;
		return $self;
	}

	sub initialize {
		my ( $self, ) = @_;
		$self->memory(15);
		$self->program->additional_params(
			[ "--variant", $self->in, "--out", $self->out ] );
		$self->output_by_type( 'vcf', $self->out );
		$self->output_by_type( 'idx', $self->idx_from_vcf( $self->out ) );
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
		$self->memory(16);
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
		$self->output_by_type( 'idx', $self->idx_from_vcf($output) );
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
		$self->memory(12);
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
		$self->output_by_type( 'vcf', $output );
		$self->output_by_type( 'idx', $self->idx_from_vcf($output) );
	}
	1;
}
#######################################################
{

	package BedToolsJob;
	use Data::Dumper;
	use Program;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self =
		  $class->SUPER::new( %params, program => new BedToolsProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'BEDTOOLS'} );
		$self->program->name( $class . "" );
		return $self;
	}

	1;
}
#######################################################
{

	package intersectBed;
	use Data::Dumper;
	use Program;
	our @ISA = qw( BedToolsJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;

		$self->program->basic_params( [ "> " . $self->out, ] );
		return $self;
	}

	1;
}
#######################################################
{

	package closestBed;
	use Data::Dumper;
	use Program;
	our @ISA = qw( BedToolsJob );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new( %params, );
		bless $self, $class;

		$self->program->basic_params( [ "> " . $self->out, ] );
		$self->output_by_type( 'vcf', $self->out );
		return $self;
	}

	1;
}

#######################################################
{

	package SamtoolsJob;
	use Program;
	our @ISA = qw( Job );

	sub new {
		my ( $class, %params ) = @_;
		my $self =
		  $class->SUPER::new( %params, program => new SamtoolsProgram() );
		bless $self, $class;
		$self->program->path( $self->project()->{'CONFIG'}->{'SAMTOOLS'} );
		$self->program->name("samtools");
		return $self;
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

	sub bai_from_bam {
		my ( $self, $bam ) = @_;
		my $bai = $bam;
		$bai =~ s/bam$/bai/;
		return $bai;
	}

	1;
}
#######################################################
{
	use Data::Dumper;

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
		my $aligned;
		if ( $lane->{PL} =~ m/illumina/i ) {
			my @sai;
			if ( $lane->{FORWARD} || $lane->{BAM1} ) {
				my $type = $lane->{BAM1} ? 'BAM1' : 'FORWARD';
				my $align = Align->new(
					params   => $params,
					previous => [$self],
					lane     => $lane,
					type     => $type,
				);
				$align->align_fastq();

				push( @sai, $align );
			}
			if ( $lane->{REVERSE} || $lane->{BAM2} ) {
				my $type = $lane->{BAM2} ? 'BAM2' : 'REVERSE';
				my $align = Align->new(
					params   => $params,
					previous => [$self],
					lane     => $lane,
					type     => $type,
				);
				$align->align_fastq();
				push( @sai, $align );
			}
			if ( $lane->{BAM} ) {
				my $type  = 'BAM';
				my $align = Align->new(
					params   => $params,
					previous => [$self],
					lane     => $lane,
					type     => $type,
				);
				$align->align_fastq();
				push( @sai, $align );
			}
			$aligned = SaiToBam->new(
				lane     => $lane,
				params   => $params,
				previous => \@sai
			);
		}
		elsif ( $lane->{PL} =~ m/solid/i ) {    #
			$aligned = Bowtie->new(
				lane     => $lane,
				params   => $params,
				previous => [$self],
			);
		}
		my $sort_sam =
		  SortSam->new( params => $params, previous => [$aligned] );
		$self->{last_job} = $sort_sam;
	}
	1;
}
#######################################################
{

	package Aligner;
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

	sub platform {
		my ( $self, ) = @_;
		return $self->lane->{PL};
	}
	1;
}
#######################################################
{

	package Align;
	our @ISA = qw( Aligner );

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
		my $qsub_param   = "-pe mpi $proc";
		$self->program->name('bwa');
		my $illumina_fastq_qualities_flag = "";
		$illumina_fastq_qualities_flag = "-I" if $self->lane->{ILLUMINA_FASTQ};
		my $bam_param = "";

		if ( $type eq 'BAM1' ) {
			$bam_param = '-b1';
		}
		elsif ( $type eq 'BAM2' ) {
			$bam_param = '-b2';
		}
		elsif ( $type eq 'BAM' ) {
			$bam_param = '-b0';
		}
		$self->program->basic_params(
			[
				'aln', '-q 5', "-t $proc", $illumina_fastq_qualities_flag,
				$colorspace, "-f $out", $genome, $bam_param, $in
			]
		);
		$self->qsub_params($qsub_param);
		$self->out($out);
		$self->output_by_type( 'fastq',  $in );
		$self->output_by_type( 'genome', $genome );
		$self->do_not_delete('genome');
		$self->do_not_delete('fastq');
	}

	1;
}
#######################################################
{

	package Bowtie;
	our @ISA = qw( Aligner );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->memory(5);
		my $lane = $self->lane();
		my $out  = $self->project()->file_prefix() . "." . $lane->{ID} . '.bam';
		$self->program->path( $self->project->{CONFIG}->{BOWTIE} );
		$self->program->name('bowtie');
		my $view =
		  $self->project()->{'CONFIG'}->{'SAMTOOLS'} . "/samtools view";

		if ( $lane->{F3_CFASTA} && $lane->{R3_CFASTA} ) {
			$self->program->basic_params(
				[
					"-f " . $lane->{F3_CFASTA} . " " . $lane->{R3_CFASTA},
					"--Q1 " . $lane->{F3_QUAL},
					"--Q2 " . $lane->{R3_QUAL},
					"| $view -bS - >",
					$out
				]
			);
		}
		$self->out($out);
		$self->output_by_type( 'bam', $out );
		return $self;
	}
	1;
}
#######################################################
{

	package Tophat;
	our @ISA = qw( Aligner );

	sub new {
		my ( $class, %params ) = @_;
		my $self = $class->SUPER::new(%params);
		bless $self, $class;
		$self->memory(10);
		my $proc = 8;
		$self->processors($proc);
		my $qsub_param = "-pe mpi $proc";
		my $lane       = $self->lane();
		my $out_dir    = $self->project->dir . '/' . $lane->{ID};
		$self->project->mkdir($out_dir);
		my $out = $out_dir . '/tophat_out/accepted_hits.bam';
		$self->program->path( $self->project->{CONFIG}->{TOPHAT} );
		$self->program->name('tophat');

		$self->program->basic_params( [ "-p $proc", "-r $PI", "", "", ] );

		$self->out($out);
		$self->output_by_type( 'bam', $out );
		return $self;
	}
	1;
}
#######################################################
{
	use Data::Dumper;

	package CollectInsertSizeMetrics;
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
		my $input = $self->first_previous->output_by_type('bam');
		my $name  = $self->project()->file_prefix() . ".is";
		my $hist  = "$name.hist";
		my $out   = "$name.txt";
		my $ref   = $self->project()->{CONFIG}->{GENOME};
		$self->program->additional_params(
			[
				"HISTOGRAM_FILE=$hist",
				"METRIC_ACCUMULATION_LEVEL=READ_GROUP",
				"INPUT=$input",
				"OUTPUT=$out",
				"STOP_AFTER=10000000",
				"REFERENCE_SEQUENCE=$ref",
			]
		);
		$self->out($out);
		$self->output_by_type( 'txt',  $out );
		$self->output_by_type( 'hist', $hist );
	}
	1;
}

#######################################################
{
	use Data::Dumper;

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
				"CREATE_INDEX=true", "SORT_ORDER=coordinate",
			]
		);
		$self->out($output);
		$self->output_by_type( 'bam', $output );
		my $bai = $self->bai_from_bam($output);
		$self->output_by_type( 'bai', $bai );
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
		$self->output_by_type( 'bai', $output );
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
		my $bai      = $self->bai_from_bam($output);
		$self->output_by_type( 'bai', $bai );
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
		$self->output_by_type( 'bam',     $output );
		my $bai = $self->bai_from_bam($output);
		$self->output_by_type( 'bai', $bai );
	}
	1;
}
#######################################################
