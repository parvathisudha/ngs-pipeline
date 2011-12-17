#!/usr/bin/perl
use strict;
use lib qw(/data/software/pipeline/);
use ConfigRun;
use Project;
use TaskScheduler;
use Time::HiRes qw ( time sleep);
use Data::Dumper;
use Getopt::Long;
use XML::Simple;
use Job;
use Jobs;
use JobManager;

####### get arguments      ###
my ( $config_file, $mode, $debug, $rerun );
GetOptions(
	'config=s' => \$config_file,
	'mode=s'   => \$mode,
	'debug'    => \$debug,
	'rerun=s'  => \$rerun,         #out|done|both
);

####### general parameters ###
my $config         = XMLin("$config_file");
my $project        = Project->new( $config, $debug );
my $task_scheduler = TaskScheduler->new( $project, $debug );
my $job_manager    = JobManager->new($debug);
my $params         = {
	config    => $config,
	rerun     => $rerun,
	scheduler => $task_scheduler,
	project   => $project,
	memory    => 1,
	manager   => $job_manager,
};

# making right folder structure
$project->make_folder_structure();

##############################
system("date") unless $debug;

##############################
my $KG           = $project->{'CONFIG'}->{'KG'};
my $hapmap       = $project->{'CONFIG'}->{'HAPMAP'};
my $omni         = $project->{'CONFIG'}->{'OMNI'};
my $dbSNP        = $project->{'CONFIG'}->{'DBSNP'};
my $indels_mills = $project->{'CONFIG'}->{'INDELS_MILLS_DEVINE'};
my $cgi          = $project->{'CONFIG'}->{'CGI'};
my $eur          = $project->{'CONFIG'}->{'EURKG'};

####### Add Jobs #############
my $root_job = RootJob->new( params => $params, previous => undef );

my @lanes_processing;

for my $lane ( @{ $project->get_lanes() } ) {
	next if $lane->{PL} eq 'virtual';
	my $process_lane = ProcessLane->new(
		params   => $params,
		previous => [$root_job],
		lane     => $lane
	);
	push( @lanes_processing, $process_lane->last_job );
}

my $join_lane_bams = MergeSamFiles->new(
	params   => $params,
	previous => [@lanes_processing],
	out      => $project->file_prefix() . ".bam",
);

my $mark_duplicates = MarkDuplicates->new(
	params   => $params,
	previous => [$join_lane_bams],
);

my $bam2cfg = Bam2cfg->new(
	params   => $params,
	previous => [$mark_duplicates],
);
$bam2cfg->do_not_delete('cfg');

my $brdMax = BreakdancerMax->new(
	params   => $params,
	previous => [$bam2cfg],
);

my @chr = $project->read_intervals();
my @snps;
my @indels;
for my $chr (@chr) {
	my $realigner_target_creator = RealignerTargetCreator->new(
		params   => $params,
		interval => $chr,
		previous => [$mark_duplicates],
	);
	my $indel_realigner = IndelRealigner->new(
		params   => $params,
		interval => $chr,
		previous => [$realigner_target_creator],
	);
	my $sort_realigned =
	  SortSam->new( params => $params, previous => [$indel_realigner] );
	my $count_covariates =
	  CountCovariates->new( params => $params, previous => [$sort_realigned] );
	my $table_recalibration = TableRecalibration->new(
		params   => $params,
		previous => [$count_covariates]
	);
	my $index_recalibrated = BuildBamIndex->new(
		params   => $params,
		previous => [$table_recalibration]
	);
	my $call_snps = UnifiedGenotyper->new(
		params         => $params,
		previous       => [$index_recalibrated],
		variation_type => "SNP"
	);    #
	my $call_indels = UnifiedGenotyper->new(
		params         => $params,
		previous       => [$index_recalibrated],
		variation_type => "INDEL"
	);    #
	push( @snps,   $call_snps );
	push( @indels, $call_indels );
}

my $combine_snps = CombineVariants->new(
	out      => $project->file_prefix() . ".SNP.vcf",
	params   => $params,
	previous => \@snps
);

my $combine_indels = CombineVariants->new(
	out      => $project->file_prefix() . ".INDEL.vcf",
	params   => $params,
	previous => \@indels
);

my $snps_variant_recalibrator = VariantRecalibrator->new(
	params            => $params,
	previous          => [$combine_snps],
	additional_params => [
"--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap",
"--resource:omni,known=false,training=true,truth=false,prior=12.0 $omni",
"--resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbSNP",
"-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ",
		"-mode SNP",
	]
);

my $snps_apply_recalibration = ApplyRecalibration->new(
	params            => $params,
	previous          => [$snps_variant_recalibrator],
	additional_params => [ "-mode SNP", ]
);

my $indels_variant_recalibrator = VariantRecalibrator->new(
	params            => $params,
	previous          => [$combine_indels],
	additional_params => [
"--resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $indels_mills",
		"-an QD -an FS -an HaplotypeScore -an ReadPosRankSum",
		"-mode INDEL",
	]
);

my $indels_apply_recalibration = ApplyRecalibration->new(
	params            => $params,
	previous          => [$indels_variant_recalibrator],
	additional_params => [ "-mode INDEL", ]
);

my $variations = CombineVariants->new(
	out      => $project->file_prefix() . ".variations.vcf",
	params   => $params,
	previous => [ $snps_apply_recalibration, $indels_apply_recalibration ]
);

my $phase_variations = ReadBackedPhasing->new(
	params   => $params,
	previous => [$variations],
	bam      => $mark_duplicates->output_by_type('bam'),
);

my $filter_low_qual = FilterLowQual->new(
	params   => $params,
	previous => [$phase_variations]
);

my $variant_annotator = VariantAnnotator->new(
	additional_params => [
		"--comp:KG,VCF $KG",
		"--comp:HapMap,VCF $hapmap",
		"--comp:OMNI,VCF $omni",
		"--comp:CGI,VCF $cgi",
		"--resource:EUR_FREQ $eur",
		"-E EUR_FREQ.AF",
		"--resource:CGI_FREQ,VCF $cgi",
		"-E CGI_FREQ.AF",
		"--resource:KG_FREQ,VCF $KG",
		"-E KG_FREQ.AF",
	],
	params   => $params,
	previous => [$filter_low_qual]
);

my $effect_prediction = SnpEff->new(
	params   => $params,
	previous => [$variant_annotator],    #
);

my $effect_annotator = VariantAnnotator->new(
	additional_params => [
		"--annotation SnpEff",
		"--snpEffFile " . $effect_prediction->output_by_type('vcf'),
	],
	params   => $params,
	previous => [ $variant_annotator, $effect_prediction ]    #
);

my $constraints_out = $effect_prediction->output_by_type('vcf') . ".constraints.bed";
my $evolution_constraints = intersectBed->new(
	additional_params => [
		"-a " . $effect_prediction->output_by_type('vcf'),
		"-b " . $project->{'CONFIG'}->{'CONSTRAINTS'},
		"-u",
		"> $constraints_out",
	],
	params   => $params,
	previous => [ $effect_prediction ]    #
);
$evolution_constraints->output_by_type('vcf', $constraints_out);

my $bgzip = Bgzip->new(
	params   => $params,
	previous => [$effect_annotator]                           #
);

my $tabix = Tabix->new(
	params   => $params,
	previous => [$bgzip]                                      #
);

#result files:
$mark_duplicates->do_not_delete('metrics');
$mark_duplicates->do_not_delete('bam');
$mark_duplicates->do_not_delete('bai');

$brdMax->do_not_delete('max');
$brdMax->do_not_delete('bed');
$brdMax->do_not_delete('fastq');

$combine_snps->do_not_delete('vcf');
$combine_snps->do_not_delete('idx');
$combine_indels->do_not_delete('vcf');
$combine_indels->do_not_delete('idx');

$snps_variant_recalibrator->do_not_delete('recal_file');
$snps_variant_recalibrator->do_not_delete('tranches_file');
$snps_variant_recalibrator->do_not_delete('rscript_file');

$snps_apply_recalibration->do_not_delete('vcf');
$snps_apply_recalibration->do_not_delete('idx');

$indels_variant_recalibrator->do_not_delete('recal_file');
$indels_variant_recalibrator->do_not_delete('tranches_file');
$indels_variant_recalibrator->do_not_delete('rscript_file');

$indels_apply_recalibration->do_not_delete('vcf');
$indels_apply_recalibration->do_not_delete('idx');

$variations->do_not_delete('vcf');
$variations->do_not_delete('idx');
$phase_variations->do_not_delete('vcf');
$phase_variations->do_not_delete('idx');
$filter_low_qual->do_not_delete('vcf');
$filter_low_qual->do_not_delete('idx');
$variant_annotator->do_not_delete('vcf');
$variant_annotator->do_not_delete('idx');
$effect_prediction->do_not_delete('vcf');
$effect_prediction->do_not_delete('idx');
$effect_annotator->do_not_delete('vcf');
$effect_annotator->do_not_delete('idx');
$bgzip->do_not_delete('gz');
$tabix->do_not_delete('tbi');
$evolution_constraints->do_not_delete('vcf');

if ($mode eq 'ALL'){
	$job_manager->start();
}

if ($mode eq 'CLEAN'){
	$job_manager->clean();
}

