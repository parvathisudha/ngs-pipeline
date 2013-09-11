use strict;
use Vcf;
use Getopt::Long;
use Data::Dumper;
####### constants ############
my $excel_string_length_limit = 100;

####### get arguments      ###
my ( $in, $out, $annotation );
GetOptions(
	'in=s'         => \$in,
	'out=s'        => \$out,
	'annotation=s' => \$annotation,
);

my $vcf = Vcf->new( file => $in );
$vcf->parse_header( silent => 1 );

my @hgmd_format = qw/HGMD.HGMDID HGMD.confidence HGMD.disease/;
my @KG_format   =
  qw/KG_FREQ.AF KG_FREQ.AFR_AF KG_FREQ.AMR_AF KG_FREQ.ASN_AF KG_FREQ.EUR_AF/;
my @vep_format = split /,/, $annotation;

# Do some simple parsing. Most thorough but slowest way how to get the data.
my @result_header = (
	'CHROM',  'POS',  'ID',       'REF',
	'ALT',    'AF',   'PROGRAM',  'SVTYPE',
	'SVLEN',  'END',  @KG_format, 'CGI_FREQ',
	'FILTER', 'CASE', 'CONTROL',  @hgmd_format,
	@vep_format,
);

my $num_vep_format = scalar @vep_format;

open OUT, ">$out" or die "Can't open $out file for writing results\n";

print OUT join( "\t", @result_header ), "\n";
while ( my $x = $vcf->next_data_hash() ) {
	my $ref         = safe_excel_string( [ $x->{'REF'} ] );
	my $alt_alleles = safe_excel_string( $x->{'ALT'} );
	my $filter      = join( ',', @{ $x->{'FILTER'} } );
	my @to_print    = (
		$x->{'CHROM'},
		$x->{'POS'},
		$x->{'ID'},
		$ref,
		$alt_alleles,
		$x->{'INFO'}->{'AF'},
		$x->{'INFO'}->{'set'},
		$x->{'INFO'}->{'SVTYPE'},
		$x->{'INFO'}->{'SVLEN'},
		$x->{'INFO'}->{'END'},
		( map { $x->{'INFO'}->{$_} } @KG_format ),
		$x->{'INFO'}->{'CGI_FREQ.AF'},
		$filter,
		$x->{'INFO'}->{'CASE.set'},
		$x->{'INFO'}->{'CONTROL.set'},
		( map { $x->{'INFO'}->{$_} } @hgmd_format ),
		( map { $x->{'INFO'}->{$_} } @vep_format),
	);
	my $string_result = join( "\t", (@to_print) );
	print OUT $string_result . "\n";
}

close OUT;

sub safe_excel_string {
	my ($array) = @_;
	my $result_string = join( ',', @{$array} );
	if ( length $result_string >= $excel_string_length_limit ) {
		my $i = 0;
		$result_string = join( ',', map { $i++; "LONG_ALLELE$i" } @{$array} );
	}
	return $result_string;
}

