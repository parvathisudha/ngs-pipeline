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
my ( $config_file, $mode, $debug );
GetOptions(
	'config=s' => \$config_file,
	'mode=s'   => \$mode,
	'debug'    => \$debug
);

####### general parameters ###
my $input   = "/data/genomes/cgi/69Genomes.reference.20111024.EUR_AF.vcf";
my $config  = XMLin("$config_file");
my $project = Project->new( $config, $debug );
my $task_scheduler = TaskScheduler->new( $project, $debug );
my $job_manager    = JobManager->new();
my $params         = {
	config    => $config,
	scheduler => $task_scheduler,
	project   => $project,
	memory    => 1,
	manager   => $job_manager,
};

# making right folder structure
$project->make_folder_structure();

##############################
my $KG           = $project->{'CONFIG'}->{'1KG'};
my $hapmap       = $project->{'CONFIG'}->{'HAPMAP'};
my $omni         = $project->{'CONFIG'}->{'OMNI'};
my $dbSNP        = $project->{'CONFIG'}->{'DBSNP'};
my $indels_mills = $project->{'CONFIG'}->{'INDELS_MILLS_DEVINE'};
my $cgi          = $project->{'CONFIG'}->{'CGI'};
my $eur          = $project->{'CONFIG'}->{'EURKG'};

####### Add Jobs #############
my $root_job = RootJob->new( params => $params, previous => undef );

my @samples = qw/
  GS12889-1100-37-ASM
  GS12892-1100-37-ASM
  GS12890-1100-37-ASM
  GS12891-1100-37-ASM
  HG00731-1100-37-ASM
  HG00732-1100-37-ASM
  GS06985-1100-37-ASM
  GS06994-1100-37-ASM
  GS07357-1100-37-ASM
  GS10851-1100-37-ASM
  GS12004-1100-37-ASM
  GS18501-1100-37-ASM
  GS18502-1100-37-ASM
  GS18504-1100-37-ASM
  GS18505-1100-37-ASM
  GS18508-1100-37-ASM
  GS18517-1100-37-ASM
  GS18526-1100-37-ASM
  GS18537-1100-37-ASM
  GS18555-1100-37-ASM
  GS18558-1100-37-ASM
  GS18940-1100-37-ASM
  GS18942-1100-37-ASM
  GS18947-1100-37-ASM
  GS18956-1100-37-ASM
  GS19017-1100-37-ASM
  GS19020-1100-37-ASM
  GS19025-1100-37-ASM
  GS19026-1100-37-ASM
  GS19129-1100-37-ASM
  GS19648-1100-37-ASM
  GS19649-1100-37-ASM
  GS19669-1100-37-ASM
  GS19670-1100-37-ASM
  GS19700-1100-37-ASM
  GS19701-1100-37-ASM
  GS19703-1100-37-ASM
  GS19704-1100-37-ASM
  GS19735-1100-37-ASM
  GS19834-1100-37-ASM
  GS20502-1100-37-ASM
  GS20509-1100-37-ASM
  GS20510-1100-37-ASM
  GS20511-1100-37-ASM
  GS20845-1100-37-ASM
  GS20846-1100-37-ASM
  GS20847-1100-37-ASM
  GS20850-1100-37-ASM
  GS21732-1100-37-ASM
  GS21733-1100-37-ASM
  GS21737-1100-37-ASM
  GS21767-1100-37-ASM
  GS19238-1100-37-ASM
  GS19239-1100-37-ASM
  /;

for my $sample (@samples) {
	my $snp = SelectVariants->new(
		params   => $params,
		previous => [$root_job],
		in       => $input,
		out      => $project->dir . "/$sample.vcf",
		additional_params => [
			"-sn $sample",
			"--excludeNonVariants",
		],		
	);
	
	my $eff = SnpEff->new(
		params   => $params,
		previous => [$snp],
	);
	my $snps_annotator = VariantAnnotator->new(
		additional_params => [
			"--annotation SnpEff",
			"--snpEffFile " . $eff->output_by_type('vcf'),
		],
		params   => $params,
		previous => [ $snp, $eff ]
	);
}

 $job_manager->start();

