use strict;
use lib qw(/data/software/pipeline/);
use ConfigRun;
use Project;
use Data::Dumper;

####### general parameters ###
my $config = ConfigRun->new( $ARGV[0] );
my $project = Project->new( $config, 0 );
my $ids = $project->get_all_ids();
my $cmd = join(",", @$ids);
system("qdel $cmd");