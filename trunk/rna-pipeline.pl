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
my $job_manager    = JobManager->new();
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

####### Add Jobs #############
my $root_job = RootJob->new( params => $params, previous => undef );
my @lanes_processing;

my $lane   = shift @{ $project->get_lanes() };
my $tophat = Tophat->new(
	params   => $params,
	previous => [$root_job],
	lane     => $lane,
);
my $sort = BuildBamIndex->new(
	params   => $params,
	previous => [$tophat],
	in       => $tophat->output_by_type('bam')
);
#my $cufflinks = Cufflinks->new(
#	params   => $params,
#	previous => [$tophat],
#);

$job_manager->start();

