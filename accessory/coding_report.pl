use strict;
use Vcf;
use Getopt::Long;
use Data::Dumper;
####### get arguments      ###
my ($in, $out);
GetOptions( 'in=s' => \$in, 
			'out=s' => \$out, 
			);

my $vcf = Vcf->new( file => $in );
$vcf->parse_header( silent => 1 );
my $header     = $vcf->format_header();
my $vep_format = $1
  if $header =~
  m/Description=\"Consequence type as predicted by VEP. Format: (.+?)\"/;
my @vep_format = split( '\|', $vep_format );
my @blank_result = map {""} @vep_format;

# Do some simple parsing. Most thorough but slowest way how to get the data.
my @result_header = (
	'CHROM',               'POS',
	'ID',                  'REF',
	'ALT',                 'AF',
	'CGI_FREQ',            'KG_FREQ',
	'EUR_FREQ',            'FILTER',
	'SNPEFF_EFFECT',       'SNPEFF_FUNCTIONAL_CLASS',
	'SNPEFF_GENE_BIOTYPE', 'SNPEFF_GENE_NAME',
	'SNPEFF_IMPACT',       'SNPEFF_TRANSCRIPT_ID',
	'SNPEFF_CODON_CHANGE', 'SNPEFF_AMINO_ACID_CHANGE',
	'SNPEFF_EXON_ID',      'SET',
	@vep_format
);

my $num_vep_format = scalar @vep_format;

open OUT, ">$out" or die "Can't open $out file for writing results\n";

print OUT join( "\t", @result_header ), "\n";
my $exon_num   = get_feature_num_by_title( \@vep_format, 'EXON' );
my $intron_num = get_feature_num_by_title( \@vep_format, 'INTRON' );
while ( my $x = $vcf->next_data_hash() ) {
	my $alt_alleles = join( ',',  @{ $x->{'ALT'} } );
	my $filter      = join( ',',  @{ $x->{'FILTER'} } );
	my @csq         = split( ",", $x->{'INFO'}->{'CSQ'} );
###INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP. Format: Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|PolyPhen|SIFT|CANONICAL|EXON|INTRON|CCDS">
	my @first_to_print = (
			$x->{'CHROM'},
			$x->{'POS'},
			$x->{'ID'},
			$x->{'REF'},
			$alt_alleles,
			$x->{'INFO'}->{'AF'},
			$x->{'INFO'}->{'CGI_FREQ.AF'},
			$x->{'INFO'}->{'KG_FREQ.AF'},
			$x->{'INFO'}->{'EUR_FREQ.AF'},
			$filter,
			$x->{'INFO'}->{'SNPEFF_EFFECT'},
			$x->{'INFO'}->{'SNPEFF_FUNCTIONAL_CLASS'},
			$x->{'INFO'}->{'SNPEFF_GENE_BIOTYPE'},
			$x->{'INFO'}->{'SNPEFF_GENE_NAME'},
			$x->{'INFO'}->{'SNPEFF_IMPACT'},
			$x->{'INFO'}->{'SNPEFF_TRANSCRIPT_ID'},
			$x->{'INFO'}->{'SNPEFF_CODON_CHANGE'},
			$x->{'INFO'}->{'SNPEFF_AMINO_ACID_CHANGE'},
			$x->{'INFO'}->{'SNPEFF_EXON_ID'},
			$x->{'INFO'}->{'SET.set'},
		);
	for my $csq (@csq) {
		my @vep_effect = split( '\|', $csq );
		$vep_effect[$exon_num]   =~ s/\//|/;
		$vep_effect[$intron_num] =~ s/\//|/;
		if (scalar @vep_effect == $num_vep_format - 1 ){
			push(@vep_effect, "");
		}
		print OUT join( "\t", (@first_to_print, @vep_effect)), "\n";
	}
	unless(scalar @csq){
		print OUT join( "\t", (@first_to_print, @blank_result)), "\n";
	}
}

close OUT;

sub get_feature_num_by_title {
	my ( $titles, $name ) = @_;
	my $i = 0;
	for my $t (@$titles) {
		last if $t eq $name;
		$i++;
	}
	return $i;
}

