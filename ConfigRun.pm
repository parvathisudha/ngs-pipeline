package ConfigRun;
use strict;
use Lane;
use Data::Dumper;

#NAME=S000005.14.March.2011
#DIR=/data/results/Denis
#GENOME=/data/genomes/bwa/hg18/hg18.fasta
#ARCHIVED=0
#BWA=/data/software/bwa-0.5.8c
#SAMTOOLS=/data/software/samtools-0.1.12a
#BEDTOOLS=/data/software/BEDTools-Version-2.10.1/bin
#VCFCODINGSNPS=/data/software/vcfCodingSnps.v1.5
#/data/CURENT/101127_SN150_0129_B205D1ABXX/      2,3,4,5
#/data/CURENT/101222_SN150_0130_B203DNABXX/      1,2,3,4

sub new {
	my ( $class,  $file )  = @_;
	my ( $params, $lanes );
	$params = {};
	eval{
		( $params, $lanes ) = _read_config("params.txt", $params);
	};
	eval{
		( $params, $lanes ) = _read_config("/data/software/pipeline/params.txt", $params);
	};	
	( $params, $lanes ) = _read_config($file, $params);

	#$params->{'DIR'} = $params->{'DIR'} . '/' . $params->{'PROJECT'};
	my $self = {
		USER          => $params->{'USER'},
		EMAIL         => $params->{'EMAIL'},
		PROJECT       => $params->{'PROJECT'},
		DIR           => $params->{'DIR'},
		GENOME        => $params->{'GENOME'},
		ARCHIVED      => $params->{'ARCHIVED'},
		BWA           => $params->{'BWA'},
		SAMTOOLS      => $params->{'SAMTOOLS'},
		BEDTOOLS      => $params->{'BEDTOOLS'},
		BEDGENOME     => $params->{'BEDGENOME'},
		VCFCODINGSNPS => $params->{'VCFCODINGSNPS'},
		MPIRUN        => $params->{'MPIRUN'},
		GATK          => $params->{'GATK'},
		GATKGENOMEBED => $params->{'GATKGENOMEBED'},
		SNPROD        => $params->{'SNPROD'},
		LANES         => $lanes,
	};
	bless $self, $class;
	return $self;
}

sub _read_config {
	my ($file, $params) = @_;
	open( CONFIG, "<$file" ) or die "Can't open configuration file: $file\n";
	my @lanes;
	while ( my $str = <CONFIG> ) {
		chomp $str;
		next if $str =~ m/$\#/;
		if ( $str =~ m/(.+?)\=(.+)/ ) {
			my ( $param, $value ) = ( $1, $2 );
			$value =~ s/\/$//;

			#options:NAME,DIR,GENOME
			$params->{$param} = $value;
		}
		elsif ( $str =~ m/(.+?)\%(.+)/ ) {
			my @args = split( /\%/, $str );
			if ( scalar @args == 2 ) {
				my $dir = $args[0];
				$dir =~ s/\/$//;
				my @nums = split /,/, $args[1];
				for my $num (@nums) {
					my $lane =
					  Lane->new( $dir, $num, $params->{'ARCHIVED'} , 0);#POTENTIAL BUG
					push( @lanes, $lane );
				}
			}
		}		
		else {
			my @args = split( /\s+/, $str );
			if ( scalar @args == 2 ) {
				my $dir = $args[0];
				$dir =~ s/\/$//;
				my @nums = split /,/, $args[1];
				for my $num (@nums) {
					my $lane =
					  Lane->new( $dir, $num, $params->{'ARCHIVED'} , 1 )
					  ;    #POTENTIAL BUG
					push( @lanes, $lane );
				}
			}
		}
		
	}
	close CONFIG;

	return ( $params, \@lanes );
}

return 1;
