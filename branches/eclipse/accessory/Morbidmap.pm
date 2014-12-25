package Morbidmap;
use strict;
use Data::Dumper;

sub new {
        my ( $class, $str) = @_;
	#17,20-lyase deficiency, isolated, 202110 (3)|CYP17A1, CYP17, P450C17|609300|10q24.32

        chomp $str;
        $str =~ s/\|/\t/g;
        my ($phenotype_name, $gene_name, $gene_id, $locus) = split ("\t", $str);
	my $phenotype_id = "";
	if($phenotype_name =~ m/(\d{6})/){
		$phenotype_id = $1;	
	}
	$gene_name =~ s/\s+//g;
	my @d = split(",", $gene_name);
	$gene_name = $d[0];
        my $self = {
		phenotype_id => $phenotype_id,
		phenotype_name => $phenotype_name,
                gene_id => $gene_id,
		gene_name => $gene_name,
		locus => $locus,
        };
        bless $self, $class;
	
        return $self;
}

1;

