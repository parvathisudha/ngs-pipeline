use strict;
use Getopt::Long;

my $excel_string_length_limit = 32759;

####### get arguments      ###
my ( $in, $out, $annotation, $vcf_data );
GetOptions(
        'in=s'              => \$in,
        'out=s'     => \$out,
        'vcf_data=s' => \$vcf_data,
        'annotation=s' => \$annotation,
);
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  S000010
#gff
#1       .       miRNA_primary_transcript        30366   30503   .       +       .       ID=MI0006363_1;accession_number=MI0006363;Name=hsa-mir-1302-2

my @vcf_pre_header = qw /CHROM POS ID REF ALT QUAL FILTER/;
my @vcf_header = split(",", $vcf_data);
my @ann_pre_header = qw /CHROM SOURCE FEATURE START END SCORE STRAND FRAME/;
my @ann_header = split(",", $annotation);

open IN, $in or die "Can't open $in file\n";
open OUT, ">$out" or die "Can't open $out file\n";
print OUT join("\t", (@vcf_pre_header, @vcf_header, @ann_pre_header, @ann_header )), "\n";
while(<IN>){
        chomp;
        my @data = split("\t");
        print OUT join( "\t", ( @data[0..2],
                safe_excel_string($data[3]),
                safe_excel_string($data[4]),
                @data[5..6],
                get_info_from_comments($data[7],\@vcf_header),
                @data[10..17],
                get_info_from_comments($data[18],\@ann_header),
                )), "\n";
}
close IN;
close OUT;

sub safe_excel_string {
        my ($string) = @_;
        my $result_string = $string;
        my @array = split/,/;
        if ( length $string >= $excel_string_length_limit ) {
                my $i = 0;
                $result_string = join( ',', map { $i++; "LONG_ALLELE$i" } @array );
        }
        return $result_string;
}

sub get_info_from_comments{
        my ($comment, $info_headers) = @_;
        my $data = {};
        $comment =~ s/\s+/=/g;
        my @fields = split(";", $comment);
        for (@fields){
                my ($tag, $value) = ($1, $2) if m/(.+)=(.+)/;
                $value =~ s/^"//;
                $value =~ s/"$//;
                $data->{$tag} = $value;
        }
        my @result = map{$data->{$_}} @$info_headers;
        return @result;
}