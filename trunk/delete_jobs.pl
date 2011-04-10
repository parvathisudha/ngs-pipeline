use strict;
use lib qw(/data/software/pipeline/);
use ConfigRun;
use Project;
use Data::Dumper;

####### general parameters ###
my $config = ConfigRun->new( $ARGV[0] );
my $project = Project->new( $config, 0 );

print Dumper $project->get_all_ids();