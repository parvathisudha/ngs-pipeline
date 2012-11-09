use strict;
use Getopt::Long;
use Data::Dumper;

my ($out, $chrLenFile, $chrLenFile,$coefficientOfVariation,$ploidy,$outputDir,$GCcontentProfile,$mateFile, $gemMappabilityFile);
GetOptions( 
            'out' => \$out,
            'chrLenFile' => \$chrLenFile,
            'gemMappabilityFile' => \$gemMappabilityFile,
            'coefficientOfVariation' => \$coefficientOfVariation,
            'ploidy' => \$ploidy,
            'outputDir' => \$outputDir,
            'GCcontentProfile' => \$GCcontentProfile,
            'mateFile' => \$mateFile,
        );
        
my $config = <<CONF;
[general]

chrLenFile = $chrLenFile
gemMappabilityFile = $gemMappabilityFile
coefficientOfVariation = $coefficientOfVariation
ploidy = $ploidy
GCcontentProfile = $GCcontentProfile
outputDir = $outputDir


[sample]
mateFile = $mateFile
inputFormat = BAM
mateOrientation = FR

CONF

open OUT, ">$out" or die "Can't ope file: $out\n";
print OUT $config;
close OUT;