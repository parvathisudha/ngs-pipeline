#!/usr/bin/perl
use strict;
use lib qw(/data/software/pipeline/);
use ConfigRun;
use Project;
use TaskScheduler;
use Time::HiRes qw ( time sleep);
use Data::Dumper;

####### general parameters ###
my $config = ConfigRun->new( $ARGV[0] );
my $debug  = $ARGV[1] ? $ARGV[1] : 0;
my $project        = Project->new( $config,        $debug );

my @to_del;
for my $lane ( @{ $config->{LANES} } ) {
	align( $project, $lane );    #tested
	#sai files	
	del($project->reverse_sai($lane));
	del($project->forward_sai($lane)) if $lane->{'PAIRED'};
	#unsorted bam file for lane
	del($project->bam($lane));
	#sorted bam for lane
	my $sorted_prefix = $project->sorted_prefix($lane);
	return 1 if ( -e "$sorted_prefix.bam" );
	
	sort_bam( $project, $lane );    #rewrite with picard
	index_bam( $project, $lane );   #rewrite with picard
}

merge_bams($project);                   #rewrite with picard
index_merged($project);                 #tested
mark_duplicates($project);              #tested
my @chr = $project->read_intervals();
for my $chr (@chr) {
	realigner_target_creator( $project, $chr );    #..
	indel_realigner( $project, $chr );             #..
	index_realigned( $project, $chr );             #..
	count_covariates( $project, $chr );            #..
	table_recalibration( $project, $chr );         #..
	index_recalibrated( $project, $chr );          #..
	parallel_call_SNPs( $project, $chr );          #..
	parallel_call_indels( $project, $chr );        #..
	indel_annotator( $project, $chr );           #..
	#	parallel_predict_effect( $project, $chr );     #..
	#	parallel_predict_indels_effect( $project, $chr );    #..
	variant_annotator( $project, $chr );           #..
	filter_snps( $project, $chr );                 #tested
}

my @chr = $project->read_intervals();
#merge SNPs after GATK filtration
my @snps_to_merge;
for my $chr (@chr) {
	push( @snps_to_merge, $project->filter_snps($chr) );
}
my @before_snps_merge;
for my $chr (@chr) {
	push( @before_snps_merge,
		$project->filter_snps_id($chr) );
}
my $merged_snps    = $project->merge_snps();
my $merge_snps_job = $project->merge_snps_id();

merge_parallel_vcf( $project, \@snps_to_merge, $sample_id, $merged_snps, $merge_snps_job,
	\@before_snps_merge );

#merge indels after calling
my @indels_to_merge;
for my $chr (@chr) {
	push( @indels_to_merge, $project->indel_annotator($chr) );
}
my @before_indels_merge;
for my $chr (@chr) {
	push( @before_indels_merge,
		 $project->indel_annotator_id($chr)  );
}
my $merged_indels    = $project->merge_indels();
my $merge_indels_job = $project->merge_indels_id();

merge_parallel_vcf( $project, \@indels_to_merge, $sample_id, $merged_indels, $merge_indels_job,
	\@before_indels_merge );



variant_recalibrator($project);    #
apply_recalibration($project);     #

my $recalibrated_snps     = $project->apply_recalibration();
my $recalibrated_snps_job = $project->apply_recalibration_id();

#merging recalibrated bams
my $merged_recalibrated_bam     = $project->file_prefix() . ".recal.bam";
my $merge_recalibrated_bams_job =
  'merge.r.' . $project->_get_id($merged_recalibrated_bam);
merge_bam_files( recalibrated_bams(), $merged_recalibrated_bam,
	$merge_recalibrated_bams_job, indexed_recalibrated_bams_job_names() );

#depth_coverage - TESTED
my $genome_coverage_file = $project->file_prefix() . ".cov";
my $genome_coverage_job  =
  'coverage.' . $project->_get_id($genome_coverage_file);
calculate_genome_coverage(
	$merged_recalibrated_bam, $genome_coverage_file,
	$genome_coverage_job, [$merge_recalibrated_bams_job]
);

#coverage_cumulative 
my $coverage_cumulative_file = $project->file_prefix() . ".cum";
my $coverage_cumulative_job  =
  'cum.' . $project->_get_id($coverage_cumulative_file);
coverage_cumulative(
	$genome_coverage_file, $coverage_cumulative_file,
	$coverage_cumulative_job, [$genome_coverage_job]
);



my $b_genome_coverage_file = $project->file_prefix() . ".bcov";
my $b_genome_coverage_job  =
  'b.coverage.' . $project->_get_id($b_genome_coverage_file);
calculate_bga_coverage(
	$merged_recalibrated_bam, $b_genome_coverage_file,
	$b_genome_coverage_job, [$merge_recalibrated_bams_job]
);

#callable_loci - TESTED
my $genome_callable_file    = $project->file_prefix() . ".callable";
my $genome_callable_summary = $project->file_prefix() . ".callable_s";
my $genome_bed              = $project->{CONFIG}->{GATKGENOMEBED};
my $genome_callable_job     =
  'callable.' . $project->_get_id($genome_callable_file);
callable_loci( $merged_recalibrated_bam, $genome_bed, $genome_callable_file,
	$genome_callable_summary, $genome_callable_job,
	[$merge_recalibrated_bams_job] );

#gatk depth_coverage - TESTED
my $gatk_cov           = $project->file_prefix() . ".gatkcov";
my $depth_coverage_job = 'gcov.' . $project->_get_id($gatk_cov);
depth_coverage( $merged_recalibrated_bam, $genome_bed, $gatk_cov,
	$depth_coverage_job, [$merge_recalibrated_bams_job] );

#merging snps with indels - TESTED

my $snps_with_indels     = $project->file_prefix() . ".variants.vcf";
my $snps_with_indels_job = 'merg.var.' . $project->_get_id($snps_with_indels);
my $before_merge         = [ $merge_indels_job, $recalibrated_snps_job ];
merge_parallel_vcf( $project, [ $recalibrated_snps, $merged_indels ],
	$sample_id, $snps_with_indels, $snps_with_indels_job, $before_merge,
	$snps_with_indels_job );

#selecting PASS variants - !!!!!!!!!!!!!!!!!!
my $snps_with_indels_pass     = $project->file_prefix() . ".variants.pass.vcf";
my $snps_with_indels_pass_job =
  'pass.var.' . $project->_get_id($snps_with_indels_pass);
filter_pass( $snps_with_indels, $sample_id, $snps_with_indels_pass,
	$snps_with_indels_pass_job, [$snps_with_indels_job] );

#variant evaluation - !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
my $var_stat           = $project->file_prefix() . ".stat";
my $var_evaluation_job = 'var_eval.' . $project->_get_id($var_stat);
variant_evaluation( $snps_with_indels_pass, $var_stat, $var_evaluation_job,
	[$snps_with_indels_pass_job] );

#effect prediction - TESTED
my $eff_html    = $project->file_prefix() . ".eff.html";
my $eff_vcf     = $project->file_prefix() . ".eff.vcf";
my $var_eff_job = 'var_eff.' . $project->_get_id($eff_vcf);
snpeff( $snps_with_indels_pass, $eff_html, $eff_vcf, $var_eff_job,
	[$snps_with_indels_pass_job] );

#zipping and indexing file with annotated variations - TESTED
my $zipped_vars      = $eff_vcf . '.gz';
my $indexed_vars     = $zipped_vars . '.tbi';
my $zipping_vars_job = 'bgzip.' . $project->_get_id($zipped_vars);
my $tabix_vars_job   = 'tabix.' . $project->_get_id($indexed_vars);
bgzip_file( $eff_vcf, $zipping_vars_job, [$var_eff_job] );
tabix_file( $zipped_vars, $tabix_vars_job, [$zipping_vars_job] );





breakdancer_cfg($project);#install on cluster
#breakdancer_mini($project);#install on cluster
breakdancer_max($project);#install on cluster

#clean($project);

#get_target_SNPs($project);
#calculate_coverage($project);
#write_report($project);

sub del{
	my ($file) = @_;
	push (@to_del, $file);
}
