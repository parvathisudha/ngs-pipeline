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
my $task_scheduler = TaskScheduler->new( $project, $debug );

#my $email = $config->{'EMAIL'};
#my $config_dir = $config->{'DIR'};
my $proc         = 4;    #number of processos for running bwa
my $sleep_time   = 0;
my $bwa_priority = 10;

# making right folder structure
my $script_dir = $config->{'DIR'} . '/tasks';
$project->make_folder_structure();
##############################

#bwa commands
my $bwa    = $config->{'BWA'} . "/bwa";
my $mpirun = $config->{'MPIRUN'} . "/mpirun -np $proc";
my $align  = "$bwa aln -t $proc";
my $sampe  = "$bwa sampe";
my $samse  = "$bwa samse";

#samtools commands
my $samtools = $config->{'SAMTOOLS'} . "/samtools";
my $import   = "$samtools import";
my $sort     = "$samtools sort";
my $index    = "$samtools index";
my $merge    = "$samtools merge";
my $view     = "$samtools view";
my $gatk     = $config->{'GATK'} . "/GenomeAnalysisTK.jar";
my $mark_dup = $project->{'CONFIG'}->{'PICARD'} . '/' . "MarkDuplicates.jar";
my $picard_sort_index = $project->{'CONFIG'}->{'PICARD'} . '/' . "SortSam.jar";
my $picard_index = $project->{'CONFIG'}->{'PICARD'} . '/' . "BuildBamIndex.jar";
my $call         = "";
my $genome       = $config->{'GENOME'};
my $gene_list    = $config->{'GENELIST'};
my $effect       =
  $config->{'VCFCODINGSNPS'} . "/vcfCodingSnps.v1.5 -r $genome -g $gene_list";
my $filter_interesting = "perl /data/software/filter_interesting.pl";
my $genome_coverage    = $config->{'BEDTOOLS'} . "/genomeCoverageBed";

my $break_dancer_dir = $project->{'CONFIG'}->{'BREAKDANCER'};
my $bam2cfg          = "perl $break_dancer_dir/bam2cfg.pl";
my $BreakDancerMax   = "perl $break_dancer_dir/BreakDancerMax.pl";
my $BreakDancerMini  = "perl $break_dancer_dir/BreakDancerMini.pl";
####### commands to execute ##

#define_done_jobs($project);

my $output_bam = $project->merged_sorted;
unless ( -e $output_bam ) {
	for my $lane ( @{ $config->{LANES} } ) {
		align( $project, $lane );    #tested
		sai_to_sam( $project, $lane );    #tested
		     #import_sam( $project, $lane );#to deletion
		sort_bam( $project, $lane );    #rewrite with picard
		index_bam( $project, $lane );   #rewrite with picard
	}
}

merge_bams($project);                   #rewrite with picard

#sort_merged($project);#to deletion
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
	parallel_predict_effect( $project, $chr );     #..
	parallel_predict_indels_effect( $project, $chr );    #..
	variant_annotator( $project, $chr );                 #..
	filter_snps( $project, $chr );                       #..
}

merge_snps($project);                                    #!!
merge_indels($project);                                  #!!

variant_recalibrator($project);                          #
apply_recalibration($project);                           #

#merge_vcf($project);#method to deletion

callable_loci($project);    #rewrite to use recalibrated bam

#depth_coverage($project);#rewrite to use recalibrated bam
calculate_genome_coverage($project);    #rewrite to use recalibrated bam
calculate_bga_coverage($project);       #rewrite to use recalibrated bam

bgzip($project);
tabix($project);
################################################
#FUNCTION TEMPLATE
#name($project, $input_file1, ..., $ouput_file2,..., $id, $after);
################################################
#sub bgzip_file {
#	my ($project, $vcf, ..., $zipped,..., $id, $after) = @_;
#	sleep($sleep_time);
#	return 1 if ( -e "$vcf.gz" );
#	my $program    = "bgzip $vcf";
#	$task_scheduler->submit( $id, "", $program, $after);
#}



#breakdancer_cfg($project);#install on cluster
#breakdancer_mini($project);#install on cluster
#breakdancer_max($project);#install on cluster

#clean($project);

#get_target_SNPs($project);
#calculate_coverage($project);
#write_report($project);

sub align {
	my ( $project, $lane ) = @_;
	sleep($sleep_time);

	submit_alignment(
		$project->reverse_name($lane),
		$project->reverse_sai($lane),
		$project->reverse_align_id($lane)
	);

	if ( $lane->{'PAIRED'} ) {
		submit_alignment(
			$project->forward_name($lane),
			$project->forward_sai($lane),
			$project->forward_align_id($lane)
		);
	}
}

sub sai_to_sam {
	my ( $project, $lane ) = @_;
	sleep($sleep_time);
	my $forward_reads = $project->forward_name($lane);
	my $reverse_reads = $project->reverse_name($lane);
	my $forward_sai   = $project->forward_sai($lane);
	my $reverse_sai   = $project->reverse_sai($lane);
	my $sam           = $project->sam($lane);
	my $bam           = $project->bam($lane);
	my $rg            = $lane->{'RG'};
	return 1 if ( -e $bam );
	my $program    = "";
	my $qsub_param = "";

	if ( $lane->{'PAIRED'} ) {
		$program =
"$sampe $genome $forward_sai $reverse_sai $forward_reads $reverse_reads -r $rg | $view -bS - > $bam";
		$qsub_param =
		    '-hold_jid '
		  . $project->task_id( $project->forward_align_id($lane) ) . ','
		  . $project->task_id( $project->reverse_align_id($lane) );
	}
	else {
		$program =
"$samse $genome $reverse_sai $reverse_reads -r $rg | $view -bS - > $bam";
		$qsub_param =
		  '-hold_jid ' . $project->task_id( $project->reverse_align_id($lane) );
	}
	$task_scheduler->submit( $project->sam_id($lane), $qsub_param, $program );
}

#sub import_sam {
#	my ( $project, $lane ) = @_;
#	sleep($sleep_time);
#	my $sam        = $project->sam($lane);
#	my $bam        = $project->bam($lane);
#	return 1 if (-e $bam);
#	my $program    = "$import $genome.fai $sam $bam";
#	my $qsub_param =
#	  '-hold_jid ' . $project->task_id( $project->sam_id($lane) );
#	$task_scheduler->submit( $project->import_id($lane), $qsub_param,
#		$program );
#}

sub sort_bam {
	my ( $project, $lane ) = @_;
	sleep($sleep_time);
	my $bam           = $project->bam($lane);
	my $sorted_prefix = $project->sorted_prefix($lane);
	return 1 if ( -e "$sorted_prefix.bam" );
	my $program    = "$sort $bam $sorted_prefix";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->sam_id($lane) );
	$task_scheduler->submit( $project->sorted_id($lane), $qsub_param,
		$program );
}

sub index_bam {
	my ( $project, $lane ) = @_;
	sleep($sleep_time);
	my $sorted = $project->sorted($lane);
	return 1 if ( -e $sorted );
	my $program    = "$index $sorted";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->sorted_id($lane) );
	$task_scheduler->submit( $project->index_id($lane), $qsub_param, $program );
}

sub submit_alignment {
	my ( $in, $out, $job_id ) = @_;
	sleep($sleep_time);
	return 1 if ( -e $out );
	my $program    = "$align -f $out $genome $in";
	my $qsub_param = "-pe mpi $proc -p $bwa_priority";
	$task_scheduler->submit( $job_id, $qsub_param, $program );
}

sub merge_bams {
	my ($project) = @_;
	sleep($sleep_time);

	my $output_bam = $project->merged_sorted;
	return 1 if ( -e $output_bam );

	my $lanes = $project->{'CONFIG'}->{'LANES'};
	my @lane_bams;
	for my $lane (@$lanes) {
		my $file = $project->sorted($lane);
		push( @lane_bams, $file );
	}
	my @input_bams = map( "INPUT=$_", @lane_bams );
	my $tmp_dir    = $project->dir;
	my $program    =
	    "java -Xmx10g -jar "
	  . $project->{'CONFIG'}->{'PICARD'}
	  . "/MergeSamFiles.jar "
	  . join( " ", @input_bams )
	  . " OUTPUT=$output_bam VALIDATION_STRINGENCY=SILENT TMP_DIR=$tmp_dir MAX_RECORDS_IN_RAM=2500000";

	my $qsub_param = '-hold_jid ' . $project->all_indexed_ids();
	$task_scheduler->submit( $project->merged_id(), $qsub_param, $program , 10);
}

#old vesrion
#sub merge_bams {
#	my ($project) = @_;
#	sleep($sleep_time);
#	my $lanes = $project->{'CONFIG'}->{'LANES'};
#	my @lane_bams;
#	for my $lane (@$lanes) {
#		my $file = $project->bam($lane);
#		push( @lane_bams, $file );
#	}
#	my $all_bams = join( ' ', @lane_bams );
#	return 1 if (-e $project->merged());
#	my $program = "$merge " . $project->merged() . " $all_bams";
#	if ( scalar @lane_bams == 1 ) {
#		my $lb = $lane_bams[0];
#		$program = "/bin/cp $lb " . $project->merged();
#	}
#	my $qsub_param = '-hold_jid ' . $project->all_indexed_ids();
#	$task_scheduler->submit( $project->merged_id(), $qsub_param, $program );
#}

#sub sort_merged {
#	my ($project) = @_;
#	sleep($sleep_time);
#	my $merged        = $project->merged();
#	my $sorted_prefix = $project->merged_sorted_prefix();
#	return 1 if (-e "$sorted_prefix.bam");
#	my $program       = "$sort $merged $sorted_prefix";
#	my $qsub_param = '-hold_jid ' . $project->task_id( $project->merged_id() );
#	$task_scheduler->submit( $project->merged_sorted_id(),
#		$qsub_param, $program );
#}

sub index_merged {
	my ($project) = @_;
	sleep($sleep_time);
	my $merged  = $project->merged_sorted();
	my $program = "$index $merged";
	return 1 if ( -e "$merged.bai" );
	my $qsub_param = '-hold_jid ' . $project->task_id( $project->merged_id() );
	$task_scheduler->submit( $project->merged_indexed_id(),
		$qsub_param, $program );
}

sub mark_duplicates {
	my ($project) = @_;
	sleep($sleep_time);
	my $marked = $project->mark_duplicates();
	my $merged = $project->merged_sorted();
	return 1 if ( -e $marked );
	my $tmp_dir = $project->dir;
	my $program = <<PROGRAM;
java -jar $mark_dup \\
INPUT=$merged \\
OUTPUT=$marked \\
METRICS_FILE=$marked.metrics \\
CREATE_INDEX=true \\
VALIDATION_STRINGENCY=SILENT \\
TMP_DIR=$tmp_dir
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->mark_duplicates_id(),
		$qsub_param, $program , 1);
}

sub realigner_target_creator {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $realigner_target_created = $project->realigner_target_creator($chr);
	return 1 if ( -e $realigner_target_created );
	my $marked  = $project->mark_duplicates();
	my $program = <<PROGRAM;
java -Xmx1g -jar $gatk \\
-T RealignerTargetCreator \\
-R $genome \\
-o $realigner_target_created \\
-I $marked \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->mark_duplicates_id() );
	$task_scheduler->submit( $project->realigner_target_creator_id($chr),
		$qsub_param, $program , 1);
}

sub indel_realigner {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $indel_realigned = $project->indel_realigner($chr);

	return 1 if ( -e $indel_realigned );
	my $marked                   = $project->mark_duplicates();
	my $realigner_target_created = $project->realigner_target_creator($chr);

	my $program = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T IndelRealigner \\
-I $marked \\
-R $genome \\
-targetIntervals $realigner_target_created \\
-o $indel_realigned \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid '
	  . $project->task_id( $project->realigner_target_creator_id($chr) );
	$task_scheduler->submit( $project->indel_realigner_id($chr),
		$qsub_param, $program , 4);
}

sub index_realigned {
	my ( $project, $chr ) = @_;
	my $sorted = $project->index_realigned($chr);
	return 1 if ( -e $sorted );
	sleep($sleep_time);
	my $indel_realigned = $project->indel_realigner($chr);
	my $tmp_dir         = $project->dir;
	my $program         = <<PROGRAM;
java -Xmx5g -jar $picard_sort_index \\
INPUT=$indel_realigned \\
OUTPUT=$sorted \\
CREATE_INDEX=true \\
SORT_ORDER=coordinate \\
TMP_DIR=$tmp_dir \\
VALIDATION_STRINGENCY=SILENT \\
MAX_RECORDS_IN_RAM=1250000
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->indel_realigner_id($chr) );
	$task_scheduler->submit( $project->index_realigned_id($chr),
		$qsub_param, $program , 5);
}

sub count_covariates {
	my ( $project, $chr ) = @_;
	my $recal_file = $project->count_covariates($chr);
	sleep($sleep_time);
	return 1 if ( -e $recal_file );
	my $indel_realigned = $project->index_realigned($chr);
	my $dbSNP           = $project->{'CONFIG'}->{'DBSNP'};

	my $program = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-l INFO \\
-T CountCovariates \\
-R $genome \\
-I $indel_realigned \\
-cov ReadGroupCovariate \\
-cov QualityScoreCovariate \\
-cov CycleCovariate \\
-cov DinucCovariate \\
-recalFile $recal_file \\
-B:dbsnp,VCF $dbSNP \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->index_realigned_id($chr) );
	$task_scheduler->submit( $project->count_covariates_id($chr),
		$qsub_param, $program , 4);
}

sub table_recalibration {
	my ( $project, $chr ) = @_;
	my $bam_recal = $project->table_recalibration($chr);
	sleep($sleep_time);
	return 1 if ( -e $bam_recal );
	my $recal_file      = $project->count_covariates($chr);
	my $indel_realigned = $project->index_realigned($chr);
	my $program         = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-l INFO \\
-T TableRecalibration \\
-R $genome \\
-I $indel_realigned \\
-recalFile $recal_file \\
-o $bam_recal \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->count_covariates_id($chr) );
	$task_scheduler->submit( $project->table_recalibration_id($chr),
		$qsub_param, $program , 4);
}

sub index_recalibrated {
	my ( $project, $chr ) = @_;
	my $index_recal = $project->index_recalibrated($chr);
	sleep($sleep_time);
	return 1 if ( -e $index_recal );
	my $bam_recal = $project->table_recalibration($chr);
	my $tmp_dir   = $project->dir;
	my $program   = <<PROGRAM;
java -jar $picard_index \\
INPUT=$bam_recal \\
OUTPUT=$index_recal \\
TMP_DIR=$tmp_dir
PROGRAM

	my $qsub_param =
	  '-hold_jid '
	  . $project->task_id( $project->table_recalibration_id($chr) );
	$task_scheduler->submit( $project->index_recalibrated_id($chr),
		$qsub_param, $program , 1);
}

sub parallel_call_SNPs {
	my ( $project, $chr ) = @_;
	my $bam_recal = $project->table_recalibration($chr);
	sleep($sleep_time);
	my $recal_file = $project->count_covariates($chr);
	my $gatk_vcf   = $project->parallel_gatk_vcf($chr);
	return 1 if ( -e $gatk_vcf );
	my $id              = $project->{'CONFIG'}->{'PROJECT'};
	my $dbSNP           = $project->{'CONFIG'}->{'DBSNP'};
	my $stand_call_conf = $project->{'CONFIG'}->{'GATK_stand_call_conf'};
	my $stand_emit_conf = $project->{'CONFIG'}->{'GATK_stand_emit_conf'};
	my $program         = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-R $genome \\
-T UnifiedGenotyper \\
-I $bam_recal \\
-B:dbsnp,VCF $dbSNP \\
-o $gatk_vcf \\
-stand_call_conf $stand_call_conf \\
-stand_emit_conf $stand_emit_conf \\
-dcov 80 -U \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->index_recalibrated_id($chr) );
	$task_scheduler->submit( $project->parallel_gatk_vcf_id($chr),
		$qsub_param, $program , 4);
}

sub parallel_call_indels {
	my ( $project, $chr ) = @_;
	my $bam_recal = $project->table_recalibration($chr);
	sleep($sleep_time);
	my $gatk_vcf = $project->parallel_call_indels($chr);
	return 1 if ( -e $gatk_vcf );
	my $id              = $project->{'CONFIG'}->{'PROJECT'};
	my $dbSNP           = $project->{'CONFIG'}->{'DBSNP'};
	my $stand_call_conf = $project->{'CONFIG'}->{'GATK_stand_call_conf'};
	my $stand_emit_conf = $project->{'CONFIG'}->{'GATK_stand_emit_conf'};
	my $program         = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-R $genome \\
-T UnifiedGenotyper \\
-glm INDEL \\
-I $bam_recal \\
-B:dbsnp,VCF $dbSNP \\
-o $gatk_vcf \\
-stand_call_conf $stand_call_conf \\
-stand_emit_conf $stand_emit_conf \\
-dcov 80 -U \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->index_recalibrated_id($chr) );
	$task_scheduler->submit( $project->parallel_call_indels_id($chr),
		$qsub_param, $program , 4);
}

sub parallel_predict_effect {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $eff_vcf  = $project->parallel_predict_effect($chr);
	my $gatk_vcf = $project->parallel_gatk_vcf($chr);
	return 1 if ( -e $eff_vcf );
	my $program    = "$effect -s $gatk_vcf -o $eff_vcf -l $eff_vcf.log";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->parallel_gatk_vcf_id($chr) );
	$task_scheduler->submit( $project->parallel_predict_effect_id($chr),
		$qsub_param, $program );
}

sub parallel_predict_indels_effect {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $eff_vcf  = $project->parallel_predict_indels_effect($chr);    #
	my $gatk_vcf = $project->parallel_call_indels($chr);
	return 1 if ( -e $eff_vcf );
	my $program    = "$effect -s $gatk_vcf -o $eff_vcf -l $eff_vcf.log";
	my $qsub_param =
	  '-hold_jid '
	  . $project->task_id( $project->parallel_call_indels_id($chr) );
	$task_scheduler->submit( $project->parallel_predict_indels_effect_id($chr),
		$qsub_param, $program );
}

sub merge_parallel_vcf {
	my ( $project, $vcfs, $merged_name, $after_id, $id ) = @_;
	my $merged_vcf = $merged_name;
	return 1 if ( -e $merged_vcf );
	sleep($sleep_time);
	my @vcfs = map { "-B:sample,VCF $_" } @$vcfs;
	my $all_vcf = join( ' ', @vcfs );
	my $program = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T CombineVariants \\
-R $genome \\
-o $merged_vcf \\
$all_vcf \\
-variantMergeOptions UNION \\
--assumeIdenticalSamples
PROGRAM
	my $qsub_param = '-hold_jid ' . $after_id;
	$task_scheduler->submit( $id, $qsub_param, $program , 4);
}

sub merge_snps {
	my ($project) = @_;
	my $merged_vcf = $project->merge_vcf;
	merge_parallel_vcf( $project, $project->all_snps_eff_vcf,
		$project->merge_vcf, $project->all_annotated, $project->merge_vcf_id );
}

sub merge_indels {
	my ($project) = @_;
	my $merged_vcf = $project->merge_indels;
	merge_parallel_vcf(
		$project, $project->all_indels_eff_vcf,
		$project->merge_indels, $project->all_indels_annotated,
		$project->merge_indels_id
	);
}

sub merge_vcf {
	my ($project) = @_;
	my $merged_vcf = $project->merge_vcf();
	return 1 if ( -e $merged_vcf );
	sleep($sleep_time);
	my @vcfs = map { "-B:sample,VCF $_" } $project->all_gatk_vcf();
	my $all_vcf = join( ' ', @vcfs );
	my $program = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T CombineVariants \\
-R $genome \\
-o $merged_vcf \\
$all_vcf \\
-variantMergeOptions UNION \\
--assumeIdenticalSamples
PROGRAM
	my $qsub_param = '-hold_jid ' . $project->all_annotated();
	$task_scheduler->submit( $project->merge_vcf_id(), $qsub_param, $program , 4);
}

sub variant_recalibrator {
	my ($project) = @_;
	sleep($sleep_time);
	my $variant_recalibrator = $project->variant_recalibrator();
	return 1 if ( -e $variant_recalibrator );
	my $dbSNP       = $project->{'CONFIG'}->{'DBSNP'};
	my $hapmap      = $project->{'CONFIG'}->{'HAPMAP'};
	my $omni        = $project->{'CONFIG'}->{'OMNI'};
	my $merged_snps = $project->merge_snps();
	my $program     = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T VariantRecalibrator  \\
-R $genome \\
-B:input,VCF $merged_snps \\
-B:hapmap,VCF,known=false,training=true,truth=true,prior=15.0 $hapmap \\
-B:omni,VCF,known=false,training=true,truth=false,prior=12.0 $omni \\
-B:dbsnp,VCF,known=true,training=false,truth=false,prior=8.0 $dbSNP \\
-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an HRun \\
-recalFile $variant_recalibrator.recal \\
-tranchesFile $variant_recalibrator \\
-rscriptFile $variant_recalibrator.plots.R
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merge_snps_id() );
	$task_scheduler->submit( $project->variant_recalibrator_id(),
		$qsub_param, $program , 4);
}

sub apply_recalibration {
	my ($project) = @_;
	sleep($sleep_time);
	my $recalibrated_vcf = $project->apply_recalibration();
	return 1 if ( -e $recalibrated_vcf );
	my $merged_snps          = $project->merge_snps();
	my $variant_recalibrator = $project->variant_recalibrator();
	my $program              = <<PROGRAM;
java -Xmx3g -jar $gatk \\
-T ApplyRecalibration  \\
-R $genome \\
-B:input,VCF $merged_snps \\
--ts_filter_level 99.0 \\
-recalFile $variant_recalibrator.recal \\
-tranchesFile $variant_recalibrator \\
-o $recalibrated_vcf
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->variant_recalibrator_id() );
	$task_scheduler->submit( $project->apply_recalibration_id(),
		$qsub_param, $program , 3);
}

sub breakdancer_cfg {
	my ($project) = @_;
	sleep($sleep_time);
	my $merged                 = $project->merged_sorted();
	my $breakdancer_cfg_result = $project->breakdancer_cfg();

	#return 1 if (-e $breakdancer_cfg_result);
	my $program    = "$bam2cfg $merged > $breakdancer_cfg_result";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->breakdancer_cfg_id(),
		$qsub_param, $program );
}

sub breakdancer_max {
	my ($project) = @_;
	sleep($sleep_time);
	my $breakdancer_cfg_result = $project->breakdancer_cfg();
	my $breakdancer_max_result = $project->breakdancer_max();

	#return 1 if (-e $breakdancer_max_result);
	my $program =
	  "$BreakDancerMax $breakdancer_cfg_result > $breakdancer_max_result";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->breakdancer_cfg_id() );
	$task_scheduler->submit( $project->breakdancer_max_id(),
		$qsub_param, $program );
}

sub breakdancer_mini {
	my ($project) = @_;
	sleep($sleep_time);
	my $breakdancer_cfg_result  = $project->breakdancer_cfg();
	my $breakdancer_mini_result = $project->breakdancer_mini();

	#return 1 if (-e $breakdancer_mini_result);
	my $program =
	  "$BreakDancerMini $breakdancer_cfg_result > $breakdancer_mini_result";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->breakdancer_cfg_id() );
	$task_scheduler->submit( $project->breakdancer_mini_id(),
		$qsub_param, $program );
}

sub define_done_jobs {
	my ($project) = @_;
	sleep($sleep_time);
	print join( ",", $project->get_all_written_files() ), "\n";
}

sub callable_loci {
	my ($project) = @_;
	sleep($sleep_time);
	my $bam_recal    = $project->table_recalibration();
	my $loci         = $project->callable_loci();
	my $loci_summary = $project->callable_loci_summary();
	my $id           = $project->{'CONFIG'}->{'PROJECT'};
	my $genome_bed   = $project->{'CONFIG'}->{'GATKGENOMEBED'};
	return 1 if ( -e $loci );
	my $program = <<PROGRAM;
java -Xmx4g -jar $gatk -R $genome -T CallableLoci -I $bam_recal -o $loci -l INFO -format BED -maxDepth 160 -L $genome_bed -summary $loci_summary
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->callable_loci_id(),
		$qsub_param, $program , 4);
}

sub depth_coverage {
	my ($project) = @_;
	sleep($sleep_time);
	my $merged = $project->merged_sorted();
	my $out    = $project->depth_coverage();
	return 1 if ( -e $out );
	my $genome_bed = $project->{'CONFIG'}->{'GATKGENOMEBED'};
	my $program    = <<PROGRAM;
java -jar $gatk -R $genome -T DepthOfCoverage -I $merged -L $genome_bed -o $out
PROGRAM
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->depth_coverage_id(),
		$qsub_param, $program , 1);
}

sub predict_effect {
	my ($project) = @_;
	sleep($sleep_time);
	my $gatk_vcf = $project->gatk_vcf();
	my $eff_vcf  = $project->eff_vcf();
	return 1 if ( -e $eff_vcf );
	my $program    = "$effect -s $gatk_vcf -o $eff_vcf -l $eff_vcf.log";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->gatk_vcf_id() );
	$task_scheduler->submit( $project->eff_vcf_id(), $qsub_param, $program );
}

sub bgzip {
	my ($project) = @_;
	sleep($sleep_time);
	my $eff_vcf = $project->apply_recalibration();
	return 1 if ( -e $project->bgzip() );
	my $program    = "bgzip $eff_vcf";
	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->apply_recalibration_id() );
	$task_scheduler->submit( $project->bgzip_id(), $qsub_param, $program );
}

sub tabix {
	my ($project) = @_;
	sleep($sleep_time);
	my $bgz = $project->bgzip();
	return 1 if ( -e $project->tabix() );
	my $program    = "tabix -p vcf $bgz";
	my $qsub_param = '-hold_jid ' . $project->task_id( $project->bgzip_id() );
	$task_scheduler->submit( $project->tabix_id(), $qsub_param, $program );
}

sub variant_annotator {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $annotated = $project->variant_annotator($chr);
	return 1 if ( -e $annotated );
	my $snps      = $project->parallel_predict_effect($chr);
	my $bam_recal = $project->table_recalibration($chr);
	my $snps_1KG  = $project->{'CONFIG'}->{'1KG'};
	my $hapmap    = $project->{'CONFIG'}->{'HAPMAP'};
	my $omni      = $project->{'CONFIG'}->{'OMNI'};
	my $dbSNP     = $project->{'CONFIG'}->{'DBSNP'};
	my $program   = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T VariantAnnotator \\
-l INFO \\
-R $genome \\
-I $bam_recal \\
-o $annotated \\
--useAllAnnotations \\
-B:variant,VCF $snps \\
-B:comp1KG,VCF $snps_1KG \\
-B:compHapMap,VCF $hapmap \\
-B:compOMNI,VCF $omni \\
-B:compDBSNP,VCF $dbSNP \\
-L $chr
PROGRAM
	my $qsub_param =
	  '-hold_jid '
	  . $project->task_id( $project->parallel_predict_effect_id($chr) );
	$task_scheduler->submit( $project->variant_annotator_id($chr),
		$qsub_param, $program , 4);

}

sub filter_snps {
	my ( $project, $chr ) = @_;
	sleep($sleep_time);
	my $filtered = $project->filter_snps($chr);
	return 1 if ( -e $filtered );
	my $annotated = $project->variant_annotator($chr);
	my $indels    = $project->parallel_predict_indels_effect($chr);
	my $program   = <<PROGRAM;
java -Xmx4g -jar $gatk \\
-T VariantFiltration \\
-R $genome \\
-o $filtered \\
-B:variant,VCF $annotated \\
-B:mask,VCF $indels \\
--maskName InDel \\
--clusterWindowSize 10 \\
--filterExpression "AB < 0.2 || MQ0 > 50" \\
--filterName "20110516" \\
-L $chr
PROGRAM
	my $after =
	    $project->task_id( $project->variant_annotator_id($chr) ) . ","
	  . $project->task_id( $project->parallel_predict_indels_effect_id($chr) );
	my $qsub_param = '-hold_jid ' . $after;
	$task_scheduler->submit( $project->filter_snps_id($chr),
		$qsub_param, $program , 4);
}

sub calculate_genome_coverage {
	my ($project) = @_;
	sleep($sleep_time);
	my $bam_recal = $project->table_recalibration();
	return 1 if ( -e $project->genome_coverage() );
	my $coverage_file = $project->genome_coverage();
	my $program       =
	    "$genome_coverage -ibam $bam_recal -g "
	  . $project->{'CONFIG'}->{'BEDGENOME'}
	  . " > $coverage_file";

	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->genome_coverage_id(),
		$qsub_param, $program );
}

sub calculate_bga_coverage {
	my ($project) = @_;
	sleep($sleep_time);
	my $bam_recal = $project->table_recalibration();
	return 1 if ( -e $project->genome_coverage_bga() );
	my $coverage_file = $project->genome_coverage_bga();

	my $program =
	    "$genome_coverage -bga -ibam $bam_recal -g "
	  . $project->{'CONFIG'}->{'BEDGENOME'}
	  . " > $coverage_file";

	my $qsub_param =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->genome_coverage_bga_id(),
		$qsub_param, $program );
}

#sub move_bedtools_results {
#	my ($project) = @_;
#	sleep($sleep_time);
#	return 1 if (-e $project->genome_coverage_bga());
#	return 1 if (-e $project->genome_coverage());
#	my $out_dir    = $project->output_dir();
#	my $genome_cov = $out_dir . '/' . $project->genome_coverage_id() . "*";
#	my $bga_cov    = $out_dir . '/' . $project->genome_coverage_bga_id() . "*";
#	my $program    =
#	    "/bin/cp $genome_cov "
#	  . $project->genome_coverage() . "\n"
#	  . "/bin/cp $bga_cov "
#	  . $project->genome_coverage_bga() . "\n";
#	my $qsub_param =
#	    '-hold_jid '
#	  . $project->task_id( $project->genome_coverage_id() ) . ','
#	  . $project->task_id( $project->genome_coverage_bga_id() );
#	$task_scheduler->submit( $project->move_bedtools_results_id(),
#		$qsub_param, $program );
#}

sub clean() {
	my ($project) = @_;
	sleep($sleep_time);
	my $files_array_to_remove = $project->get_garbage_files();
	my $files_to_remove       = join( ' ', @$files_array_to_remove );
	my $program               = "/bin/rm -f $files_to_remove";
	my $qsub_param            =
	  '-hold_jid ' . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->clean_id(), $qsub_param, $program );
}

sub call_SNPs {
	my ($project) = @_;
	sleep($sleep_time);
	my $merged   = $project->merged_sorted();
	my $gatk_vcf = $project->gatk_vcf();
	my $id       = $project->{'CONFIG'}->{'PROJECT'};
	my $dbSNP    = $project->{'CONFIG'}->{'DBSNP'};
	return 1 if ( -e $gatk_vcf );
	my $nump       = 5;
	my $memory     = $nump * 4;
	my $memory_str = '-Xmx' . $memory . 'g';
	my $memory_qs  = 'mem_total=' . $memory . 'G';
	my $program    = <<PROGRAM;
java $memory_str -jar $gatk -R $genome -T UnifiedGenotyper -I $merged -B:dbsnp,VCF $dbSNP -o $gatk_vcf -nt $nump \\
-stand_call_conf 50.0 \\
-stand_emit_conf 10.0 \\
-dcov 80 -U
PROGRAM
	my $qsub_param =
	    "-pe mpi $nump -l $memory_qs "
	  . '-hold_jid '
	  . $project->task_id( $project->merged_indexed_id() );
	$task_scheduler->submit( $project->gatk_vcf_id(), $qsub_param, $program );
}

#
#sub get_target_SNPs{
#	my ($project) = @_;
#	my $sorted = $project->sorted($lane);
#	my $program = "$index $sorted";
#	my $qsub_param = '-hold_jid ' . $project->sorted_id($lane);
#	$task_scheduler->submit($config->{'USER'}, $project->index_id($lane), $qsub_param, $program, $email);
#
#}

#
#sub write_report{
#	my ($project) = @_;
#	my $sorted = $project->sorted($lane);
#	my $program = "$index $sorted";
#	my $qsub_param = '-hold_jid ' . $project->sorted_id($lane);
#	$task_scheduler->submit($config->{'USER'}, $project->index_id($lane), $qsub_param, $program, $email);
#
#}
