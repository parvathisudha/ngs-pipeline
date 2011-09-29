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
my $config         = XMLin("$config_file");
my $project        = Project->new( $config, $debug );
my $task_scheduler = TaskScheduler->new( $project, $debug );
my $job_manager    = JobManager->new();
my $params         = {
	config    => $config,
	scheduler => $task_scheduler,
	project   => $project,
	memory    => 1,
	manager => $job_manager,
};

# making right folder structure
$project->make_folder_structure();
####### Add Jobs #############
my $root_job = RootJob->new( params => $params, previous => undef);
my @lanes_processing;
for my $lane ( @{ $project->get_lanes() } ) {
	my $process_lane = ProcessLane->new( params => $params, previous => [$root_job] , lane => $lane );
	push( @lanes_processing, $process_lane );
}

$job_manager->start();
