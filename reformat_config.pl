#use strict;
#open IN, $ARGV[0] or die "Can't open file $ARGV[0]\n";
#
#while (my $str = <IN>) {
#	chomp $str;
#	unless($str =~ m/\t/ || $str =~ m/\%/){
#		print "$str\n";
#		next;
#	}
#	if ( $str =~ m/(.+?)\%(.+)/ ) {
#		my @args = split( /\%/, $str );
#		if ( scalar @args == 3 ) {
#			my $dir = $args[0];
#			$dir =~ s/\/$//;
#			my @nums = split /,/, $args[1];
#			my $rg = $args[2];
#			for my $num (@nums) {
#				my $lane =
#				  Lane->new( $dir, $num, $params->{'ARCHIVED'}, 0, $rg )
#				  ;    #POTENTIAL BUG
#				push( @lanes, $lane );
#			}
#		}
#	}
#	else {
#		my @args = split( /\t+/, $str );
#		if ( scalar @args == 3 ) {
#			my $dir = $args[0];
#			$dir =~ s/\/$//;
#			my @nums = split /,/, $args[1];
#			my $rg = $args[2];
#			for my $num (@nums) {
#				my $lane =
#				  Lane->new( $dir, $num, $params->{'ARCHIVED'}, 1, $rg )
#				  ;    #POTENTIAL BUG
#				push( @lanes, $lane );
#			}
#		}
#	}
#
#}
#close IN;
