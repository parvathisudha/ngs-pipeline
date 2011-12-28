use strict;

my ($in, $out, $max_kg_freq, $max_eur_freq, $max_cgi_freq) = @ARGV;

open IN, "<$in";
open OUT, ">$out";
while(<IN>){
        if(m/^#/){
                print OUT;
        }
        else{
                my $kg = $1 if m/KG_FREQ\.AF=(.+?)\;/;
                my $eur = $1 if m/EUR_FREQ\.AF=(.+?)\;/;
                my $cgi = $1 if m/CGI_FREQ\.AF=(.+?)\;/;
                if($kg < $max_kg_freq && $eur < $max_eur_freq && $cgi < $max_cgi_freq){
                        print OUT;
                }
        }
}
close IN;
close OUT;