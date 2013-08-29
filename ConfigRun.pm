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
	my $self = $params;
	$self->{'LANES'} = $lanes;
	bless $self, $class;
	return $self;
}

return 1;
