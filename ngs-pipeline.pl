#!/usr/bin/perl
use strict;
use lib qw(/data/software/pipeline/);
use ConfigRun;
use Project;
use TaskScheduler;
use Time::HiRes qw ( time sleep);
use Data::Dumper;
use Getopt::Long;
use XML::Simple;
use Job;
use Jobs;
use JobManager;

####### get arguments      ###
my ( $config_file, $mode, $debug, $rerun );
Getopt::Long::Configure("pass_through");
GetOptions(
	'config=s' => \$config_file,
	'mode=s'   => \$mode,
	'debug'    => \$debug,
	'rerun=s'  => \$rerun,         #out|done|both
);

if ( scalar @ARGV ) {
	warn "Unknown options:", Dumper \@ARGV, "\n";
	exit 0;
}

####### general parameters ###
my $params_file = "params.xml";
my $config = XMLin( "$config_file", ForceArray => [ 'LANE', ], );
$0 =~ /^(.+[\\\/])[^\\\/]+[\\\/]*$/;
my $path = $1 || "./";
$path =~ s/\/$//;
my $params = XMLin("$path/$params_file");
for my $param ( keys %{ $params->{PARAMETERS} } ) {
	$config->{PARAMETERS}->{$param} = $params->{PARAMETERS}->{$param}
	  unless exists $config->{PARAMETERS}->{$param};
}

my $project        = Project->new( $config,        $debug );
my $task_scheduler = TaskScheduler->new( $project, $debug );
my $job_manager    = JobManager->new($debug);
my $params         = {
	config    => $config,
	rerun     => $rerun,
	scheduler => $task_scheduler,
	project   => $project,
	memory    => 1,
	manager   => $job_manager,
};

# making right folder structure
$project->make_folder_structure();

##############################
system("date") unless $debug;

##############################
my $KG                = $project->{'CONFIG'}->{'KG'};
my $SVKG              = $project->{'CONFIG'}->{'SVKG'};
my $HGMD              = $project->{'CONFIG'}->{'HGMD'};
my $hapmap            = $project->{'CONFIG'}->{'HAPMAP'};
my $omni              = $project->{'CONFIG'}->{'OMNI'};
my $dbSNP             = $project->{'CONFIG'}->{'DBSNP'};
my $indels_mills      = $project->{'CONFIG'}->{'INDELS_MILLS_DEVINE'};
my $cgi               = $project->{'CONFIG'}->{'CGI'};
my $eur               = $project->{'CONFIG'}->{'EURKG'};
my $group_vcf         = $project->{'CONFIG'}->{'CASE'};
my $control_group_vcf = $project->{'CONFIG'}->{'CONTROL'};
my $max_freq          = $project->{'CONFIG'}->{'MAXFREQ'};
my $loci              = $project->{'CONFIG'}->{'LOCI'};

#############################

my @snpeff_coding_classes = qw/
  MODERATE
  HIGH
  SYNONYMOUS_START
  NON_SYNONYMOUS_START
  START_GAINED
  NON_SYNONYMOUS_STOP
  /;

my @vep_coding_classes = qw/
  FRAMESHIFT_CODING
  NON_SYNONYMOUS_CODING
  NON_SYNONYMOUS_CODING
  SPLICE_SITE
  STOP_GAINED
  STOP_LOST
  SPLICE_SITE
  ESSENTIAL_SPLICE_SITE
  COMPLEX_IN
  PARTIAL_CODON
  /;
my @vep_severe = qw/
	transcript_ablation
	splice_donor_variant
	splice_acceptor_variant
	stop_gained
	frameshift_variant
	stop_lost
	initiator_codon_variant
	inframe_insertion
	inframe_deletion
	missense_variant
	transcript_amplification
	splice_region_variant
	incomplete_terminal_codon_variant
	mature_miRNA_variant
	TFBS_ablation
	TFBS_amplification
	TF_binding_site_variant
	regulatory_region_variant
	regulatory_region_ablation
	regulatory_region_amplification
  /; 
my @vep_non_severe = qw/
	frameshift_variant&non_coding_exon_variant
	frameshift_variant&3_prime_UTR_variant
	frameshift_variant&5_prime_UTR_variant
	frameshift_variant&NMD_transcript_variant
	NMD_transcript_variant
	nc_transcript_variant
  /;     
my $coding_classes_string        = join( '|', @vep_coding_classes );
my $vep_severe_string        = join( '|', @vep_severe );
my $vep_non_severe_string        = join( '|', @vep_non_severe );
my $snpeff_coding_classes_string = join( '|', @snpeff_coding_classes );
####### Add Jobs #############
my $root_job = RootJob->new( params => $params, previous => undef );

my @lanes_processing;

for my $lane ( @{ $project->get_lanes() } ) {
	next if $lane->{PL} eq 'virtual';
	my $process_lane = ProcessLane->new(
		params   => $params,
		previous => [$root_job],
		lane     => $lane
	);
	my $telomere_count = CalcTelomeres->new(
		params   => $params,
		previous => [$root_job],
		lane     => $lane
	);
	my $telomere_norm = NormalizeTelomeres->new(
		params   => $params,
		previous => [$telomere_count],
	);
	push( @lanes_processing, $process_lane->last_job );
}

my $join_lane_bams = MergeSamFiles->new(
	params   => $params,
	previous => [@lanes_processing],
	out      => $project->file_prefix() . ".bam",
);

my $mark_duplicates = MarkDuplicates->new(
	params   => $params,
	previous => [$join_lane_bams],
);

my $insert_size_metrics = CollectInsertSizeMetrics->new(
	params   => $params,
	previous => [$mark_duplicates],
);
$insert_size_metrics->do_not_delete('txt');
$insert_size_metrics->do_not_delete('hist');

my $coverage = DepthOfCoverage->new(
	params            => $params,
	previous          => [$mark_duplicates],
	additional_params => [
		"--omitDepthOutputAtEachBase", "--omitIntervalStatistics",
		"--omitLocusTable",
	]
);
$coverage->do_not_delete('sample_summary');
$coverage->do_not_delete('sample_statistics');
$coverage->do_not_delete('sample_summary');

#------------- BreakDancer -----------------------
my $bam2cfg;
my $brdMax;
if ($project->{'CONFIG'}->{'EXPERIMENT'} eq 'WGS'){
	$bam2cfg = Bam2cfg->new(
		params   => $params,
		previous => [$mark_duplicates],
	);
	$bam2cfg->do_not_delete('cfg');
	
	$brdMax = BreakdancerMax->new(
		params   => $params,
		previous => [$bam2cfg],
	);	
	$brdMax->do_not_delete('max');
	$brdMax->do_not_delete('bed');
	$brdMax->do_not_delete('fastq');
}

#------------- Copy Number Variations ------------

runCNVSeq( $mark_duplicates, $params );
runFREEC( $mark_duplicates, $params );

#------------- Pindel ----------------------------
my $dedup_index_link = Ln->new(
	params   => $params,
	previous => [$mark_duplicates],
	in       => $mark_duplicates->output_by_type('bai'),
	out      => $mark_duplicates->output_by_type('bam') . ".bai",
);

$dedup_index_link->do_not_delete('link');

my $pindel_config = PindelConfig->new(
	params            => $params,
	previous          => [$insert_size_metrics],
	additional_params => [
		"--bam",             $mark_duplicates->output_by_type('bam'),
		"--id",              $project->sample_id,
		"--min_insert_size", $project->{'CONFIG'}->{'PINDEL_MIN_INSERT_SIZE'},
	]
);
$pindel_config->do_not_delete('cfg');
my $pindel;

if($project->{'CONFIG'}->{'EXPERIMENT'} eq 'WGS'){
	my $pindel = Pindel->new(
		params            => $params,
		previous          => [$pindel_config],
		additional_params => [ "--breakdancer", $brdMax->out, ]
	);	
}
else{
	my $pindel = Pindel->new(
		params            => $params,
		previous          => [$pindel_config],
	);	
}


$pindel->do_not_delete('SV_D');
$pindel->do_not_delete('SV_INV');
$pindel->do_not_delete('SV_LI');
$pindel->do_not_delete('SV_SI');
$pindel->do_not_delete('SV_TD');
$pindel->do_not_delete('SV_BP');

#------------- Process Pindel output ----------------
my $pindel_results = {};
for my $pindel_out ( @{ $pindel->VEP_compatible_files } ) {
	my $pindel2vcf = Pindel2Vcf->new(
		params   => $params,
		previous => [$pindel],
		in       => $pindel_out,
		out      => $pindel_out . ".vcf",
	);
	$pindel2vcf->do_not_delete('vcf');
	my $pindel_left_aligned = LeftAlignVariants->new(
		params   => $params,
		previous => [$pindel2vcf],
		in       => $pindel2vcf->out,
		out      => $pindel2vcf->out . ".la.vcf",
	);
	$pindel_left_aligned->do_not_delete('vcf');
	$pindel_left_aligned->do_not_delete('idx');

	#------------- ChipSeq analysis
	if ( $project->{'CONFIG'}->{'CHIPSEQ'} ) {
		my $chipseq = IntersectVcfBed->new(
			out => $pindel_left_aligned->output_by_type('vcf') . ".chipseq.vcf",
			bed => $project->{'CONFIG'}->{'CHIPSEQ'},
			params   => $params,
			previous => [$pindel_left_aligned]    #
		);
		$chipseq->do_not_delete('vcf');
		$chipseq->do_not_delete('idx');

		my $chipseq_idx = IGVTools->new(
			out               => $chipseq->output_by_type('vcf') . ".idx",
			params            => $params,
			previous          => [$chipseq],
			additional_params => [ "index", $chipseq->output_by_type('vcf'), ]
		);
		$chipseq_idx->do_not_delete('main');
	}

	$pindel_results->{$pindel_out} = $pindel_left_aligned;
}

my @pindel_for_VEP =
  map { $pindel_results->{$_} } @{ $pindel->VEP_compatible_files };
my @pindel_for_SnpEff =
  map { $pindel_results->{$_} } @{ $pindel->SnpEff_compatible_files };

#------------- Merge Pindel output ----------------
my $combined_pindel = CombineVariants->new(
	out      => $project->file_prefix() . ".PINDEL.vcf",
	params   => $params,
	previous => \@pindel_for_VEP,
);
$combined_pindel->do_not_delete('vcf');
$combined_pindel->do_not_delete('idx');

#------------- GATK SNP and INDEL calling --------

my @chr = $project->read_intervals();
my @snps;
my @indels;
for my $chr (@chr) {
	my $realigner_target_creator = RealignerTargetCreator->new(
		params   => $params,
		interval => $chr,
		previous => [$mark_duplicates],
	);
	my $indel_realigner = IndelRealigner->new(
		params   => $params,
		previous => [$realigner_target_creator],
		interval => $chr,
	);
	my $sort_realigned =
	  SortSam->new( params => $params, previous => [$indel_realigner] );
	my $count_covariates =
	  CountCovariates->new( params => $params, previous => [$sort_realigned] );
	my $table_recalibration = TableRecalibration->new(
		params   => $params,
		previous => [$count_covariates],
		interval => $chr,
	);
	my $index_recalibrated = BuildBamIndex->new(
		params   => $params,
		previous => [$table_recalibration]
	);
	my $call_snps = UnifiedGenotyper->new(
		params         => $params,
		previous       => [$index_recalibrated],
		variation_type => "SNP",
		interval => $chr,
	);    #
	my $call_indels = UnifiedGenotyper->new(
		params         => $params,
		previous       => [$index_recalibrated],
		variation_type => "INDEL",
		interval => $chr,
	);    #
	push( @snps,   $call_snps );
	push( @indels, $call_indels );
}

my $combine_snps = CombineVariants->new(
	out      => $project->file_prefix() . ".SNP.vcf",
	params   => $params,
	previous => \@snps
);
my $combine_indels = CombineVariants->new(
	out      => $project->file_prefix() . ".INDEL.vcf",
	params   => $params,
	previous => \@indels
);

my $snps_variant_recalibrator = VariantRecalibrator->new(
	params            => $params,
	previous          => [$combine_snps],
	additional_params => [
"--resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hapmap",
"--resource:omni,known=false,training=true,truth=false,prior=12.0 $omni",
"--resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $dbSNP",
"-an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ",
		"-mode SNP",
	]
);
$snps_variant_recalibrator->do_not_delete('recal_file');
$snps_variant_recalibrator->do_not_delete('tranches_file');
$snps_variant_recalibrator->do_not_delete('rscript_file');
my $snps_apply_recalibration = ApplyRecalibration->new(
	params            => $params,
	previous          => [$snps_variant_recalibrator],
	additional_params => [ "-mode SNP", ]
);
$snps_apply_recalibration->do_not_delete('vcf');
$snps_apply_recalibration->do_not_delete('idx');
my $indels_variant_recalibrator = VariantRecalibrator->new(
	params            => $params,
	previous          => [$combine_indels],
	additional_params => [
"--resource:mills,VCF,known=true,training=true,truth=true,prior=12.0 $indels_mills",
		"-an QD -an FS -an HaplotypeScore -an ReadPosRankSum",
		"-mode INDEL",
	]
);
$indels_variant_recalibrator->do_not_delete('recal_file');
$indels_variant_recalibrator->do_not_delete('tranches_file');
$indels_variant_recalibrator->do_not_delete('rscript_file');
my $indels_apply_recalibration = ApplyRecalibration->new(
	params            => $params,
	previous          => [$indels_variant_recalibrator],
	additional_params => [ "-mode INDEL", ]
);
my $phase_variations = ReadBackedPhasing->new(
	params   => $params,
	previous => [$snps_apply_recalibration],
	bam      => $mark_duplicates->output_by_type('bam'),
);

my $variations = CombineVariants->new(
	out      => $project->file_prefix() . ".variations.vcf",
	params   => $params,
	previous => [ $phase_variations, $indels_apply_recalibration ]
);

my $filter_low_qual = FilterLowQual->new(
	params   => $params,
	previous => [$variations]
);

#Merge GATK and Pindel outputs

my $gatk_and_pindel_combined = CombineVariants->new(
	out      => $project->file_prefix() . ".gp.vcf",
	params   => $params,
	previous => [ $filter_low_qual, $combined_pindel ]
);
$gatk_and_pindel_combined->program->clean_additional_params;
$gatk_and_pindel_combined->program->additional_params(
	[
		"--variant:gatk",
		$filter_low_qual->output_by_type('vcf'),
		"--variant:pindel",
		$combined_pindel->output_by_type('vcf'),
		"-o",
		$gatk_and_pindel_combined->output_by_type('vcf'),
		"-genotypeMergeOptions PRIORITIZE",
		"-priority gatk,pindel",
	]
);
$gatk_and_pindel_combined->do_not_delete('vcf');
$gatk_and_pindel_combined->do_not_delete('idx');

my $variant_annotator = VariantAnnotator->new(
	additional_params => [
		"--comp:KG,VCF $KG",
		"--resource:KG_FREQ,VCF $KG",
		"-E KG_FREQ.AF",
		"-E KG_FREQ.AFR_AF",
		"-E KG_FREQ.AMR_AF",
		"-E KG_FREQ.ASN_AF",
		"-E KG_FREQ.EUR_AF",

		"--comp:CGI,VCF $cgi",
		"--resource:CGI_FREQ,VCF $cgi",
		"-E CGI_FREQ.AF",

		"--resource:CASE,VCF $group_vcf",
		"-E CASE.set",
		"--resource:CONTROL,VCF $control_group_vcf",
		"-E CONTROL.set",
		"--resource:HGMD,VCF $HGMD",
		"-E HGMD.HGMDID",
		"-E HGMD.confidence",
		"-E HGMD.disease",
		"-E HGMD.hyperlink",
		"-E HGMD.mutationType",
		"-E HGMD.nucleotideChange",
		"-E HGMD.pmid",
		"-E HGMD.variantType",
	],
	params   => $params,
	previous => [$gatk_and_pindel_combined]
);

##################### HGMD VARIATIONS ##############

my $hgmd_vcf_grep = GrepVcf->new(
	params       => $params,
	basic_params => ["--regexp HGMDID"],
	previous     => [$variant_annotator]    #
);
$hgmd_vcf_grep->do_not_delete('vcf');
$hgmd_vcf_grep->do_not_delete('idx');

my $hgmd_vcf_table = VariantsToTable->new(
	params            => $params,
	out               => $project->file_prefix() . ".hgmd.txt",
	additional_params => [
		"-F CHROM -F POS -F ID -F REF -F ALT -F AF -F CGI_FREQ\.AF",
		"-F KG_FREQ\.AF -F EUR_FREQ\.AF -F QUAL -F FILTER",
"-F HGMD\.HGMDID -F HGMD\.confidence -F HGMD\.disease -F HGMD\.hyperlink",
"-F HGMD\.mutationType -F HGMD\.nucleotideChange -F HGMD\.pmid -F HGMD\.variantType",
		"-F SNPEFF_EFFECT -F SNPEFF_FUNCTIONAL_CLASS",
		"-F SNPEFF_GENE_BIOTYPE -F SNPEFF_GENE_NAME -F SNPEFF_IMPACT",
		"-F SNPEFF_TRANSCRIPT_ID -F SNPEFF_CODON_CHANGE",
"-F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_EXON_ID -F CASE.set -F CONTROL.set",
		"--showFiltered"
	],
	previous => [$hgmd_vcf_grep]    #
);
$hgmd_vcf_table->do_not_delete('txt');
##################### ANALYSIS OF SEVERE SNVs ##############
my $rare = FilterFreq->new(
	params       => $params,
	basic_params => [ $max_freq, $max_freq, ],
	previous     => [$variant_annotator]         #
);
$rare->do_not_delete('vcf');
$rare->do_not_delete('idx');

my $rare_ann_eff = VEP->new(
	params            => $params,
	previous          => [$rare],
	out               => $rare->output_by_type('vcf') . ".vep.vcf",
	additional_params => [
		"--sift=b", "--polyphen=b", "--per_gene",
		"--hgnc",   "--ccds",       "--canonical", "--numbers",
	],
);
$rare_ann_eff->do_not_delete('vcf');
$rare_ann_eff->do_not_delete('idx');

my $coding = GrepVcf->new(
	params       => $params,
	basic_params => [ "--regexp '" . $vep_severe_string . "'" ],
	previous     => [$rare_ann_eff]                                    #
);
$coding->do_not_delete('vcf');
$coding->do_not_delete('idx');

my $rare_cod_table = CodingReport->new(
	params   => $params,
	out      => $project->file_prefix() . ".cod.txt",
	previous => [$coding]                                              #
);
$rare_cod_table->do_not_delete('txt');

my $grep_rare_cod_table = GrepTxt->new(
	params       => $params,
	basic_params => [ 
	"--regexp '" . $vep_severe_string . "'",
	"--regexp_v '" . $vep_non_severe_string . "'",
						
	 ],
	previous     => [$rare_cod_table]                                              #
);
$grep_rare_cod_table->do_not_delete('txt');

my $cod_annotate_proteins = AnnotateProteins->new(
	params            => $params,
	out               => $grep_rare_cod_table->out . '.uniprot.txt',
	additional_params => [
		"--in",
		$grep_rare_cod_table->out,
		"--id_column gene",
		"--gene_to_protein",
		$project->{'CONFIG'}->{'ENSEMBL_TO_UNIPROT'},
		"--id_type gene",
		"--uniprot_db",
		$project->{'CONFIG'}->{'UNIPROT'},
	],
	previous => [$grep_rare_cod_table]    #
);
$cod_annotate_proteins->do_not_delete('txt');

my $cod_annotate_proteins_mark = JoinTabular->new(
	params            => $params,
	out               => $cod_annotate_proteins->out . '.marked.txt',
	additional_params => [
		"--table",
		$cod_annotate_proteins->out,
		"--annotation",
		$project->{'CONFIG'}->{'GOI'},
		"--table_id_columns 22 --annotation_id_columns 0",
		"--annotation_columns 1,2,3",
		"--all_table",
	],
	previous => [ $cod_annotate_proteins, ]    #
);
$cod_annotate_proteins_mark->do_not_delete('txt');

my $loci_cod_table = AddLoci->new(
	params            => $params,
	previous          => [$cod_annotate_proteins],
	additional_params => [ "--loci $loci", ],
);
$loci_cod_table->do_not_delete('txt');

########## miRNA genes and targets   ################
my $rare_miRNA = VEP->new(
	params            => $params,
	previous          => [$rare],
	out               => $rare->output_by_type('vcf') . ".miRNA.vcf",
	additional_params => [
		"--no_consequence",
		"-custom "
		  . $project->{'CONFIG'}->{'MIRNA_SITE'}
		  . ",MIRNA_SITE,bed,overlap,0",
		"-custom "
		  . $project->{'CONFIG'}->{'MIRNA_TRANSCRIPT'}
		  . ",MIRNA_TRANSCRIPT,bed,overlap,0",
		"-custom "
		  . $project->{'CONFIG'}->{'MIRNA_MATURE'}
		  . ",MIRNA_MATURE,bed,overlap,0",
		"-custom "
		  . $project->{'CONFIG'}->{'GENES_BED'}
		  . ",GENE,bed,overlap,0",
	],
);
$rare_miRNA->do_not_delete('vcf');
$rare_miRNA->do_not_delete('idx');

### miRNA genes

my $grep_miRNA = GrepVcf->new(
	params       => $params,
	out          => $project->file_prefix() . ".mir_genes.vcf",
	basic_params => [
		"--regexp_v '" . join( '|', ( "SVTYPE=INV", "SVTYPE=RPL" ) ) . "'",
		"--regexp '"
		  . join( '|', ( "MIRNA_TRANSCRIPT", "MIRNA_MATURE" ) ) . "'",
	],
	previous => [$rare_miRNA]
);
$grep_miRNA->do_not_delete('vcf');

my $rare_miRNA_genes_report = VcfToReport->new(
	params            => $params,
	out               => $project->file_prefix() . ".mir_genes.xls",
	previous          => [$grep_miRNA],
	additional_params => [ "--annotation", "MIRNA_TRANSCRIPT,MIRNA_MATURE", ],
);
$rare_miRNA_genes_report->do_not_delete('xls');

### miRNA targets

my $grep_miRNA_targets = GrepVcf->new(
	params       => $params,
	out          => $project->file_prefix() . ".mir_targets.vcf",
	basic_params => [
		"--regexp_v '" . join( '|', ( "SVTYPE=INV", "SVTYPE=RPL" ) ) . "'",
		"--regexp '" . join( '|', ("MIRNA_SITE") ) . "'",
	],
	previous => [$rare_miRNA]
);
$grep_miRNA_targets->do_not_delete('vcf');

my $rare_miRNA_targets_report = VcfToReport->new(
	params            => $params,
	out               => $project->file_prefix() . ".mir_targets.xls",
	previous          => [$grep_miRNA_targets],
	additional_params => [ "--annotation", "MIRNA_SITE,GENE", ],
);
$rare_miRNA_targets_report->do_not_delete('xls');

my $rare_miRNA_targets_report_proteins = AnnotateProteins->new(
	params            => $params,
	out               => $project->file_prefix() . ".mir_targets.a.xls",
	additional_params => [
		"--in",
		$rare_miRNA_targets_report->out,
		"--id_column gene",
		"--gene_to_protein",
		$project->{'CONFIG'}->{'ENSEMBL_TO_UNIPROT'},
		"--id_type gene",
		"--uniprot_db",
		$project->{'CONFIG'}->{'UNIPROT'},
	],
	previous => [$rare_cod_table]    #
);
$rare_miRNA_targets_report_proteins->do_not_delete('main');

########## REGULATION ANALYSIS ######################
my $evolution_constraints_for_reg = IntersectVcfBed->new(
	out      => $variant_annotator->output_by_type('vcf') . ".constraints.vcf",
	bed      => $project->{'CONFIG'}->{'CONSTRAINTS'},
	params   => $params,
	previous => [$variant_annotator]                                           #
);
$evolution_constraints_for_reg->do_not_delete('vcf');
$evolution_constraints_for_reg->do_not_delete('idx');

my $reg_constraints_rare = FilterFreq->new(
	params       => $params,
	basic_params => [ $max_freq, $max_freq, $max_freq, ],
	previous     => [$evolution_constraints_for_reg]                           #
);
$reg_constraints_rare->do_not_delete('vcf');
$reg_constraints_rare->do_not_delete('idx');

#------------------------Fix site predictions!!!!!!!!!!!!!
my $reg_constraints_rare_table = VariantsToTable->new(
	params            => $params,
	out               => $project->file_prefix() . ".constrained.rare.txt",
	additional_params => [
		"-F CHROM -F POS -F ID -F REF -F ALT -F AF -F CGI_FREQ\.AF",
		"-F KG_FREQ\.AF -F EUR_FREQ\.AF -F QUAL",
		"-F FILTER -F SNPEFF_EFFECT -F SNPEFF_FUNCTIONAL_CLASS",
		"-F SNPEFF_GENE_BIOTYPE -F SNPEFF_GENE_NAME -F SNPEFF_IMPACT",
		"-F SNPEFF_TRANSCRIPT_ID -F SNPEFF_CODON_CHANGE",
		"-F SNPEFF_AMINO_ACID_CHANGE -F SNPEFF_EXON_ID -F CASE.set",
		"--showFiltered"
	],
	previous => [$reg_constraints_rare]    #
);
$reg_constraints_rare_table->do_not_delete('txt');

#-------------------------------------------------------------------------------------

my $in_ensemble_regulatory = GrepVcf->new(
	params       => $params,
	basic_params =>
	  [ "--regexp REGULATION", "--regexp_v '" . $coding_classes_string . "'" ],
	previous => [$reg_constraints_rare]    #
);
$in_ensemble_regulatory->do_not_delete('vcf');
$in_ensemble_regulatory->do_not_delete('idx');

my $regulatory_group_annotator = VariantAnnotator->new(
	additional_params => [
		"--resource:CASE,VCF $group_vcf",
		"-E CASE.set",
		"--resource:CONTROL,VCF $control_group_vcf",
		"-E CONTROL.set",
	],
	params   => $params,
	previous => [$in_ensemble_regulatory]
);
$regulatory_group_annotator->do_not_delete('vcf');
$regulatory_group_annotator->do_not_delete('idx');

my $near_genes = closestBed->new(
	params => $params,
	out    => $regulatory_group_annotator->output_by_type('vcf') . ".genes",
	basic_params => [
		"-t first",                                         "-a",
		$regulatory_group_annotator->output_by_type('vcf'), "-b",
		$project->{'CONFIG'}->{'TWOKBTSS'}
	],
	previous => [$regulatory_group_annotator]    #
);
$near_genes->do_not_delete('vcf');
$near_genes->do_not_delete('idx');

my $regulatory_rare_table = VariantsToTable->new(
	params            => $params,
	additional_params => [
		"-F CHROM -F POS -F ID -F REF -F ALT -F AF -F CGI_FREQ\.AF",
		"-F KG_FREQ\.AF -F EUR_FREQ\.AF -F QUAL",
		"-F FILTER -F EFF -F CASE.set -F CONTROL.set",
		"--showFiltered"
	],
	previous => [$regulatory_group_annotator]    #
);
$regulatory_rare_table->do_not_delete('txt');

my $regulatory_rare_table_with_genes = JoinTabular->new(
	params            => $params,
	out               => $regulatory_rare_table->out . '.with_genes.txt',
	additional_params => [
		"--table",      $regulatory_rare_table->out,
		"--annotation", $near_genes->out,
		"--table_id_columns 0,1,3,4 --annotation_id_columns 0,1,3,4",
		"--annotation_columns 14",
		"--annotation_header GENE_ID",
		"--table_columns 0,1,2,3,4,5,6,7,8,9,10,11,12,13",
		"--skip_annotation_header",
		"--annotation_header GENE_ID",

	],
	previous => [ $in_ensemble_regulatory, $regulatory_rare_table ]    #
);
$regulatory_rare_table_with_genes->do_not_delete('txt');

my $regulatory_rare_table_with_genes_proteins = AnnotateProteins->new(
	params => $params,
	out    => $regulatory_rare_table_with_genes->out . '.uniprot.txt',
	additional_params => [
		"--in",
		$regulatory_rare_table_with_genes->out,
		"--id_column gene",
		"--gene_to_protein",
		$project->{'CONFIG'}->{'ENSEMBL_TO_UNIPROT'},
		"--id_type gene",
		"--uniprot_db",
		$project->{'CONFIG'}->{'UNIPROT'},
	],
	previous => [$regulatory_rare_table_with_genes]    #
);
$regulatory_rare_table_with_genes_proteins->do_not_delete('txt');

my $reformat_regulation = ReformatRegulation->new(
	params            => $params,
	out               => $project->file_prefix() . ".reg.txt",
	additional_params =>
	  [ "--in", $regulatory_rare_table_with_genes->out, "--eff_column 11", ],
	previous => [$regulatory_rare_table_with_genes]    #
);
$reformat_regulation->do_not_delete('txt');

my $regulation_with_genes_marked = JoinTabular->new(
	params            => $params,
	out               => $reformat_regulation->out . '.marked.txt',
	additional_params => [
		"--table",
		$reformat_regulation->out,
		"--annotation",
		$project->{'CONFIG'}->{'GOI'},
		"--table_id_columns 14 --annotation_id_columns 0",
		"--annotation_columns 1,2,3",
		"--all_table",
	],
	previous => [ $reformat_regulation, ]    #
);
$regulation_with_genes_marked->do_not_delete('txt');

######################################################

my $bgzip = Bgzip->new(
	params   => $params,
	previous => [$variant_annotator]         #
);

my $tabix = Tabix->new(
	params   => $params,
	previous => [$bgzip]                     #
);

#result files:
$mark_duplicates->do_not_delete('metrics');
$mark_duplicates->do_not_delete('bam');
$mark_duplicates->do_not_delete('bai');



$combine_snps->do_not_delete('vcf');
$combine_snps->do_not_delete('idx');
$combine_indels->do_not_delete('vcf');
$combine_indels->do_not_delete('idx');

$indels_apply_recalibration->do_not_delete('vcf');
$indels_apply_recalibration->do_not_delete('idx');

$variations->do_not_delete('vcf');
$variations->do_not_delete('idx');
$phase_variations->do_not_delete('vcf');
$phase_variations->do_not_delete('idx');
$filter_low_qual->do_not_delete('vcf');
$filter_low_qual->do_not_delete('idx');
$variant_annotator->do_not_delete('vcf');
$variant_annotator->do_not_delete('idx');
$bgzip->do_not_delete('gz');
$tabix->do_not_delete('tbi');

#$evolution_constraints->do_not_delete('vcf');
#$effect_annotator_rare->do_not_delete('vcf');
#$constraints_rare->do_not_delete('vcf');
#$effect_annotator_rare->do_not_delete('vcf');

if ( $mode eq 'ALL' || $mode eq 'PINDEL_TRUE' ) {
	$job_manager->start();
}

if ( $mode eq 'CLEAN' ) {
	$job_manager->clean();
}

sub runCNVSeq {
	my ( $previous, $params ) = @_;
	my $hits          = $project->file_prefix() . ".hits";
	my $generate_hits = SamtoolsJob->new(
		params            => $params,
		out               => $hits,
		additional_params => [
			"view -F 4",
			$previous->output_by_type("bam"),
			"| perl -lane \'print \"\$F[2]\t\$F[3]\"\'",
			"> ", $hits,
		],
		previous => [ $previous, ]    #
	);
	$generate_hits->output_by_type( "hits", $hits );
	$generate_hits->do_not_delete('hits');

	my $cnv_res = $project->file_prefix() . ".cnv.txt";
	my $cnvseq  = CNVSeq->new(
		params            => $params,
		out               => $cnv_res,
		additional_params => [
			"--test",                 $generate_hits->output_by_type("hits"),
			"--ref",                  $project->{'CONFIG'}->{'CNVSEQ_REF'},
			"--genome",               $project->{'CONFIG'}->{'CNVSEQ_GENOME'},
			"--global-normalization", "> ",
			$cnv_res,
		],
		previous => [ $generate_hits, ]    #
	);
	$cnvseq->do_not_delete('txt');
	return $cnvseq;
}

sub runFREEC {
	my ( $previous, $params ) = @_;
	my $conf       = $project->file_prefix() . ".freec.conf";
	my $write_conf = WriteFreecConfig->new(
		params            => $params,
		out               => $conf,
		additional_params => [
			'--out'        => $conf,
			'--chrLenFile' => $project->{'CONFIG'}->{'FREEC'}->{'chrLenFile'},
			'--gemMappabilityFile' =>
			  $project->{'CONFIG'}->{'FREEC'}->{'gemMappabilityFile'},
			'--coefficientOfVariation' =>
			  $project->{'CONFIG'}->{'FREEC'}->{'coefficientOfVariation'},
			'--ploidy'    => $project->{'CONFIG'}->{'FREEC'}->{'ploidy'},
			'--outputDir' => $project->dir(),
			'--chrFiles'  => $project->{'CONFIG'}->{'FREEC'}->{'chrFiles'},
			'--mateFile'  => $previous->output_by_type('bam'),
		],
		previous => [ $previous, ]    #
	);
	$write_conf->do_not_delete('txt');

	my $freec = FREEC->new(
		params   => $params,
		previous => [ $write_conf, ]    #
	);
	$freec->do_not_delete('cnv');

	my $cnv_bed_res = $freec->out . ".bed";
	my $cnv_bed     = Cut->new(
		params            => $params,
		out               => $cnv_bed_res,
		additional_params => [ "-f 1,2,3,5", $freec->out, "> $cnv_bed_res", ],
		previous          => [ $freec, ]                                       #
	);
	$cnv_bed->do_not_delete('main');

	my $gain_out = $cnv_bed->out . ".gain.bed";
	my $gain     = Grep->new(
		params            => $params,
		out               => $gain_out,
		additional_params => [ "gain", $cnv_bed->out, "> $gain_out", ],
		previous          => [ $cnv_bed, ]                                     #
	);
	$gain->do_not_delete('main');

	my $loss_out = $cnv_bed->out . ".loss.bed";
	my $loss     = Grep->new(
		params            => $params,
		out               => $gain_out,
		additional_params => [ "loss", $cnv_bed->out, "> $loss_out", ],
		previous          => [ $cnv_bed, ]                                     #
	);
	$loss->do_not_delete('main');

	my $cnv_img_file = $freec->output_by_type('ratio') . ".png";
	my $cnv_img      = Cat->new(
		params            => $params,
		out               => $cnv_img_file,
		additional_params => [
			$project->{'CONFIG'}->{'FREEC'}->{'SCRIPTS'} . "/makeGraph.R",
			"|",
			"R --slave --args",
			$project->{'CONFIG'}->{'FREEC'}->{'ploidy'},
			$freec->output_by_type('ratio')
		],
		previous => [ $freec, ]    #
	);
	$cnv_img->do_not_delete('main');

	my $freec2circos = FREEC2Circos->new(
		params   => $params,
		previous => [ $freec, ],
	);
	$freec2circos->do_not_delete('main');

	return $freec;
}

