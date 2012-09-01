use strict;
use Data::Dumper;
use Getopt::Long;
####### get arguments      ###
my ( @file, $out, $reads_limit, $samtools, $distribution_file , $fastq_quality_filter );
GetOptions(
	'file=s@'        => \@file,
	'out=s'          => \$out,
	'reads_limit=s'  => \$reads_limit,
	'samtools=s'     => \$samtools,
#	'result=s'       => \$result_file,
	'distribution=s' => \$distribution_file,
	'fastq_quality_filter=s' => \$fastq_quality_filter,
);

my $telomere_repeats_num  = 0;
my $total_reads           = 0;
my $telomere_reads        = 0;
my $is_seq                = 0;
my $reads_len             = 0;
my $telomers_len          = 0;
my $min_tel_rep           = 10;
my $telomere_distribution = {};
my $quality_filter = "$fastq_quality_filter -q 30 -p 90";

for my $file (@file) {
	my $fh                = get_stream($file);
	my $strings_from_file = 0;
	while ( my $str = <$fh> ) {
		chomp $str;
		$strings_from_file++;
		if((($strings_from_file-2) % 4) == 0){
			my $observed_f = 0;
			my $observed_r = 0;
			$observed_f++ while ( $str =~ m/TTAGGG/g );
			$observed_r++ while ( $str =~ m/CCCTAA/g );
			$total_reads++;
			$reads_len += length $str;
			my $rep = ( $observed_f > $observed_r ) ? $observed_f : $observed_r;
			if ( exists $telomere_distribution->{$rep} ) {
				$telomere_distribution->{$rep} = $telomere_distribution->{$rep} + 1;
			}
			else {
				$telomere_distribution->{$rep} = 1;
			}
	
			if ( $rep >= $min_tel_rep ) {
				$telomers_len += length $str;
				$telomere_reads++;
				$telomere_repeats_num += $rep;
			}			
		}
		if ( $strings_from_file / 4 >= $reads_limit / ( scalar @file ) ) {
			last;
		}
	}
	close $fh;
}
#my $genome_size = 3195751584;
##http://www.ncbi.nlm.nih.gov/projects/genome/assembly/grc/human/data/index.shtml
#
#my $coverage           = $reads_len / $genome_size;
#my $tel_per_chr_length = $genome_size * $telomers_len / ( $reads_len * 46 );
#
#my $average_repeat_num_per_telomeric_region =
#  $telomere_repeats_num / $telomers_len;
#
#my $telomere_per_chr_num =
#  $average_repeat_num_per_telomeric_region * $tel_per_chr_length;
#
##my $tel_per_chr_length = $telomers_len / ($coverage * 46);
#
#open RESULT, ">$result_file" or die "Can't write to $result_file\n";
#print RESULT
#"TOTAL_READS\tTELOMERE_REPEATS\tSEQUENCED_LENGTH\tTELOMERE_SEQUENCED_LENGTH\tCOVERAGE\tTELOMERE_PER_CHR_NUM\tTELOMERE_PER_CHR_LENGTH\tTELOMERE_REP_AVERAGE\n";
#print RESULT
#"$total_reads\t$telomere_repeats_num\t$reads_len\t$telomers_len\t$coverage\t$telomere_per_chr_num\t$tel_per_chr_length\t$average_repeat_num_per_telomeric_region\n";
#close RESULT;

open DISTRIBUTION, ">$distribution_file"
  or die "Can't write to $distribution_file\n";
print DISTRIBUTION "REPEAT_LENGTH\tREADS\n" if (scalar keys %$telomere_distribution);
for ( sort { $a <=> $b } keys %$telomere_distribution ) {
	print DISTRIBUTION $_, "\t", $telomere_distribution->{$_}, "\n";
}
close DISTRIBUTION;

sub get_stream {
	my ($file) = @_;
	if ( $file =~ m/gz$/ ) {
		open FILE, "zcat $file | $quality_filter |" or die "Can't open $file\n";
		return \*FILE;
	}
	elsif ( $file =~ m/bam$/ ) {
		open FILE,
"$samtools view $file | awk \'{print \"@\"\$1\"\\n\"\$10\"\\n+\"\$1\"\\n\"\$11}\' | $quality_filter -Q 33 |"
		  or die "Can't open $file\n";
		return \*FILE;
	}
	else {
		die "Unknown file type: $file\n";
	}
}
