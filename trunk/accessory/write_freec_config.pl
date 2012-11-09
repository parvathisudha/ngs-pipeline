#!/usr/bin/perl
use strict;
use Getopt::Long;
use Data::Dumper;

my ($out, $chrLenFile, $chrLenFile,$coefficientOfVariation,$ploidy,$outputDir,$chrFiles,$mateFile, $gemMappabilityFile);
GetOptions( 
            'out=s' => \$out,
            'chrLenFile=s' => \$chrLenFile,
            'gemMappabilityFile=s' => \$gemMappabilityFile,
            'coefficientOfVariation=s' => \$coefficientOfVariation,
            'ploidy=s' => \$ploidy,
            'outputDir=s' => \$outputDir,
            'chrFiles=s' => \$chrFiles,
            'mateFile=s' => \$mateFile,
        );
        
my $config = <<CONF;
[general]

chrLenFile = $chrLenFile
gemMappabilityFile = $gemMappabilityFile
coefficientOfVariation = $coefficientOfVariation
ploidy = $ploidy
chrFiles = $chrFiles
outputDir = $outputDir


[sample]
mateFile = $mateFile
inputFormat = BAM
mateOrientation = FR

CONF

open OUT, ">$out" or die "Can't ope file: $out\n";
print OUT $config;
close OUT;