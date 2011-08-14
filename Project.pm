package Project;
use strict;
use Data::Dumper;

#file names for all tasks

sub new {
	my ( $class, $config, $debug ) = @_;
	my $time = get_time($config->{'DIR'});
	if ( $config->{'TIME'} ) {
		$time = $config->{'TIME'};
	}
	my $self = {
		CONFIG => $config,
		DEBUG  => $debug,
		TIME   => $time,
	};
	bless $self, $class;
	return $self;
}

#folders

sub dir {
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'};
}

sub script_dir {
	my ($self) = @_;
	return $self->dir() . '/tasks';
}

sub output_dir {
	my ($self) = @_;
	return $self->dir() . '/out';
}

sub error_dir {
	my ($self) = @_;
	return $self->dir() . '/errors';
}

sub ids_dir {
	my ($self) = @_;
	return $self->dir() . '/ids';
}

sub task_id_file {
	my ( $self, $job_name ) = @_;
	return $self->ids_dir() . "/$job_name.id";
}

sub task_id {
	my ( $self, $job_name ) = @_;
	my $file = $self->task_id_file($job_name);
	return $self->get_sge_id($file);
}

sub make_folder_structure {
	my ($self) = @_;
	$self->mkdir( $self->dir() );
	$self->mkdir( $self->script_dir() );
	$self->mkdir( $self->output_dir() );
	$self->mkdir( $self->error_dir() );
	$self->mkdir( $self->ids_dir() );
}

sub mkdir {
	my ( $self, $directory ) = @_;
	$self->command("mkdir $directory");
}

sub forward_name {
	my ( $self, $lane ) = @_;
	return $lane->forward_name();
}

sub reverse_name {
	my ( $self, $lane ) = @_;
	return $lane->reverse_name();
}

sub forward_prefix {
	my ( $self, $lane ) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->forward_align_id($lane);
}

sub reverse_prefix {
	my ( $self, $lane ) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->reverse_align_id($lane);
}

sub forward_sai {
	my ( $self, $lane ) = @_;
	return $self->forward_prefix($lane) . ".sai";
}

sub reverse_sai {
	my ( $self, $lane ) = @_;
	return $self->reverse_prefix($lane) . ".sai";
}

sub get_read_id {
	my ( $self, $lane, $is_forward ) = @_;
	my $read_direction = $is_forward ? 1 : 2;

	#s_$num_$read_direction_sequence.txt
	my $name =
	    $self->_get_id( $lane->{DIR} ) . '_' . 's_'
	  . $lane->{NUM} . '_'
	  . $read_direction
	  . '_sequence.txt';
	return $name;
}

sub forward_align_id {
	my ( $self, $lane ) = @_;
	my $id = 'aln.' . $self->_get_id( $self->get_read_id( $lane, 1 ) );
	return $id;
}

sub reverse_align_id {
	my ( $self, $lane ) = @_;
	return 'aln.' . $self->_get_id( $self->get_read_id( $lane, 0 ) );
}

sub sam {
	my ( $self, $lane ) = @_;
	return $self->forward_prefix($lane) . ".sam";
}

sub bam {
	my ( $self, $lane ) = @_;
	return $self->forward_prefix($lane) . ".bam";
}

sub sorted_prefix {
	my ( $self, $lane ) = @_;
	return $self->forward_prefix($lane) . ".sorted";
}

sub sorted {
	my ( $self, $lane ) = @_;
	return $self->sorted_prefix($lane) . ".bam";
}

sub file_prefix {
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->{'CONFIG'}->{'PROJECT'};
}

sub merged {
	my ($self) = @_;
	return $self->file_prefix() . ".bam";
}

sub merged_id {
	my ($self) = @_;
	return 'merge.' . $self->_get_id( $self->merged() );
}

sub merged_indexed {
	my ($self) = @_;
	return $self->file_prefix() . ".bai";
}

sub merged_indexed_id {
	my ($self) = @_;
	return 'index_merged.' . $self->_get_id( $self->merged_indexed() );
}

sub merged_sorted_prefix {
	my ($self) = @_;
	return $self->file_prefix() . ".sorted";
}

sub merged_sorted {
	my ($self) = @_;
	return $self->merged_sorted_prefix() . ".bam";
}

sub merged_sorted_id {
	my ($self) = @_;
	return 'sort_merged.' . $self->_get_id( $self->merged_sorted() );
}

sub mark_duplicates {
	my ($self) = @_;
	return $self->file_prefix() . ".dedup.bam";
}

sub mark_duplicates_id {
	my ($self) = @_;
	return 'mark_dup.' . $self->_get_id( $self->mark_duplicates() );
}

sub merge_snps {
	my ($self) = @_;
	return $self->file_prefix() . ".snps.merged.vcf";
}

sub merge_snps_id {
	my ($self) = @_;
	return 'merge.snps.' . $self->_get_id( $self->merge_snps() );
}

sub merge_vcf {
	my ($self) = @_;
	return $self->merge_snps;
}

sub merge_vcf_id {
	my ($self) = @_;
	return $self->merge_snps_id;
}

sub merge_indels {
	my ($self) = @_;
	return $self->file_prefix() . ".ind.merged.vcf";
}

sub merge_indels_id {
	my ($self) = @_;
	return 'merge.ind.' . $self->_get_id( $self->merge_indels() );
}
sub realigner_target_creator {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.intervals";
}

sub realigner_target_creator_id {
	my ($self, $chr) = @_;
	return "retar.$chr." . $self->_get_id( $self->realigner_target_creator($chr) );
}
sub indel_realigner {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.realigned.bam";
}

sub indel_realigner_id {
	my ($self,$chr) = @_;
	return "indre.$chr." . $self->_get_id( $self->indel_realigner($chr) );
}


sub index_realigned {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.realigned.sorted.bam";
}

sub index_realigned_id {
	my ($self,$chr) = @_;
	return "idxrea.$chr." . $self->_get_id( $self->index_realigned($chr) );
}

sub count_covariates {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.covar.csv";
}

sub count_covariates_id {
	my ($self, $chr) = @_;
	return "ccv.$chr." . $self->_get_id( $self->count_covariates($chr) );
}

sub table_recalibration {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.recal.bam";
}

sub table_recalibration_id {
	my ($self, $chr) = @_;
	return "tbre.$chr." . $self->_get_id( $self->table_recalibration($chr) );
}

sub index_recalibrated {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.recal.bam.bai";
}

sub index_recalibrated_id {
	my ($self, $chr) = @_;
	return "idxrec.$chr." . $self->_get_id( $self->index_recalibrated($chr) );
}

sub variant_recalibrator {
	my ($self) = @_;
	return $self->file_prefix() . ".tranches";
}

sub variant_recalibrator_id {
	my ($self, $chr) = @_;
	return "vre." . $self->_get_id( $self->variant_recalibrator() );
}

sub apply_recalibration {
	my ($self) = @_;
	return $self->file_prefix() . ".recal.vcf";
}

sub apply_recalibration_id {
	my ($self) = @_;
	return "avre." . $self->_get_id( $self->apply_recalibration() );
}

sub gatk_vcf {
	my ($self) = @_;
	return $self->file_prefix() . ".vcf";
}

sub gatk_vcf_id {
	my ($self) = @_;
	return 'gatk_snps.' . $self->_get_id( $self->gatk_vcf() );
}

sub parallel_gatk_vcf {
	my ($self,$chr) = @_;
	return $self->file_prefix() . ".$chr.vcf";
}

sub parallel_gatk_vcf_id {
	my ($self,$chr) = @_;
	return "snps.$chr." . $self->_get_id( $self->parallel_gatk_vcf($chr) );
}

sub parallel_call_indels {
	my ($self,$chr) = @_;
	return $self->file_prefix() . ".$chr.indel.vcf";
}

sub parallel_call_indels_id {
	my ($self,$chr) = @_;
	return "ind.$chr." . $self->_get_id( $self->parallel_call_indels($chr) );
}

sub parallel_predict_effect {
	my ($self,$chr) = @_;
	return $self->file_prefix() . ".$chr.eff.snp.vcf";
}

sub parallel_predict_effect_id {
	my ($self,$chr) = @_;
	return "eff.snp.$chr." . $self->_get_id( $self->parallel_predict_effect($chr) );
}

sub depth_coverage {
	my ($self) = @_;
	return $self->file_prefix() . ".depth_coverage";
}

sub depth_coverage_id {
	my ($self) = @_;
	return 'depth_cov.' . $self->_get_id( $self->depth_coverage() );
}

sub callable_loci {
	my ($self) = @_;
	return $self->file_prefix() . ".callable_loci";
}

sub callable_loci_summary {
	my ($self) = @_;
	return $self->file_prefix() . ".callable_loci_summary";
}

sub callable_loci_id {
	my ($self) = @_;
	return 'callable.' . $self->_get_id( $self->callable_loci() );
}

sub eff_vcf {
	my ($self) = @_;
	return $self->file_prefix() . ".eff.vcf";
}

sub eff_vcf_id {
	my ($self) = @_;
	return 'effect.' . $self->_get_id( $self->eff_vcf() );
}

sub bgzip {
	my ($self) = @_;
	return $self->apply_recalibration() . ".gz";
}

sub bgzip_id {
	my ($self) = @_;
	return 'bgzip.' . $self->_get_id( $self->bgzip() );
}

sub tabix {
	my ($self) = @_;
	return $self->file_prefix() . ".eff.vcf.gz.tbi";
}

sub tabix_id {
	my ($self) = @_;
	return 'tabix.' . $self->_get_id( $self->tabix() );
}

sub tmp_dir{
	my ($self) = @_;
	return $self->dir;
}

sub breakdancer_cfg {
	my ($self) = @_;
	return $self->file_prefix() . ".cfg";
}

sub breakdancer_cfg_id {
	my ($self) = @_;
	return 'brd_cfg' . $self->_get_id( $self->breakdancer_cfg() );
}

sub breakdancer_max {
	my ($self) = @_;
	return $self->file_prefix() . ".breakdancer_max";
}

sub breakdancer_max_id {
	my ($self) = @_;
	return 'brd_max' . $self->_get_id( $self->breakdancer_max() );
}

sub breakdancer_mini {
	my ($self) = @_;
	return $self->file_prefix() . ".breakdancer_mini";
}

sub breakdancer_mini_id {
	my ($self) = @_;
	return 'brd_mini' . $self->_get_id( $self->breakdancer_mini() );
}

sub variant_annotator {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.snps.annot.vcf";
}

sub variant_annotator_id {
	my ($self, $chr) = @_;
	return "a.snp.$chr." . $self->_get_id( $self->variant_annotator($chr) );
}

sub indel_annotator {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.indels.annot.vcf";
}

sub indel_annotator_id {
	my ($self, $chr) = @_;
	return "a.ind.$chr." . $self->_get_id( $self->indel_annotator($chr) );
}

sub filter_snps {
	my ($self, $chr) = @_;
	return $self->file_prefix() . ".$chr.snps.annot.filt.vcf";
}

sub filter_snps_id {
	my ($self, $chr) = @_;
	return "f.snp.$chr." . $self->_get_id( $self->filter_snps($chr) );
}

sub genome_coverage {
	my ($self) = @_;
	return $self->file_prefix() . ".genome_coverage.bed";
}

sub genome_coverage_id {
	my ($self) = @_;
	return 'g_cov.' . $self->_get_id( $self->genome_coverage() );
}

sub genome_coverage_bga {
	my ($self) = @_;
	return $self->file_prefix() . ".genome_coverage_bga.bed";
}

sub genome_coverage_bga_id {
	my ($self) = @_;
	return 'g_cov_bga.' . $self->_get_id( $self->genome_coverage_bga() );
}

sub move_bedtools_results_id {
	my ($self) = @_;
	return 'mv_bed_res.' . $self->_get_id('move_bedtools_results_id');
}

sub all_indexed_ids {
	my ($self) = @_;
	my $lanes  = $self->{'CONFIG'}->{'LANES'};
	my @ids    = "";
	for my $lane (@$lanes) {
		push( @ids, $self->task_id( $self->index_id($lane) ) );
	}
	return join( ',', @ids );
}

sub all_annotated {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->task_id( $self->filter_snps_id($chr) ) );
	}
	return join( ',', @ids );
}

sub all_filtered {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->task_id( $self->filter_snps_id($chr) ) );
	}
	return join( ',', @ids );
}

sub all_indels_annotated {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->task_id( $self->parallel_predict_indels_effect_id($chr) ) );
	}
	return join( ',', @ids );
}

sub all_gatk_vcf {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->parallel_gatk_vcf($chr) );
	}
	return @ids;
}

sub all_filtered_snps {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->filter_snps($chr) );
	}
	return \@ids;
}

sub all_snps_eff_vcf {
	my ($self) = @_;
	my @chr = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->filter_snps($chr) );
	}
	return \@ids;
}

sub all_indels_eff_vcf {
	my ($self) = @_;
	my @chr  = $self->read_intervals();
	my @ids;
	for my $chr (@chr) {
		push( @ids, $self->parallel_predict_indels_effect($chr) );
	}
	return \@ids;
}

sub sorted_id {
	my ( $self, $lane ) = @_;
	return 'sort.' . $self->_get_id( $self->sorted($lane) );
}

sub sam_id {
	my ( $self, $lane ) = @_;
	return 'sam.' . $self->_get_id( $self->sam($lane) );
}

sub import_id {
	my ( $self, $lane ) = @_;
	return 'import.' . $self->_get_id( $self->bam($lane) );
}

sub bai {
	my ( $self, $lane ) = @_;
	return $lane->forward_name() . ".sorted.bam.bai";
}

sub index_id {
	my ( $self, $lane ) = @_;
	return 'index.' . $self->_get_id( $self->bai($lane) );
}



sub get_all_written_files {
	my ($self) = @_;
	my $dir = $self->{'CONFIG'}->{'DIR'};
	my @array = <$dir/*>;
	return @array;
}

sub get_all_eff_vcf{
	my ($self) = @_;
	my @chr = $self->read_intervals();
	my @eff_files;
	for my $chr(@chr){
		push(@eff_files, $self->parallel_gatk_vcf_id($chr));
	} 
	return @chr;
}

sub get_garbage_files {
	my ($self) = @_;
	my @array;
	my $lanes = $self->{'CONFIG'}->{'LANES'};
	for my $lane (@$lanes) {
		push( @array, $self->forward_name($lane) );
		push( @array, $self->reverse_name($lane) );
		push( @array, $self->sam($lane) );
		push( @array, $self->bam($lane) );
		push( @array, $self->sorted($lane) );
		push( @array, $self->bai($lane) );

	}
	push( @array, $self->merged() );
	return \@array;
}

sub clean_id {
	my ($self) = @_;
	return 'clean.' . $self->_get_id('clean_id');
}

sub _get_id {
	my ( $self, $file ) = @_;
	$file =~ s/\//_/g;
	$file =~ s/\://g;
	$file .= "_" . $self->{'TIME'};
	return $file;
}

sub get_sge_id {
	my ( $self, $file ) = @_;
	if ( $self->{'DEBUG'} ) {
		$file = "id.txt";
	}
	open IN, "<$file" or die "Can't open file: $file for reading task id!\n";
	while (<IN>) {
		return $1 if m/Your job\s(\d+)\s/;
	}
	close IN;
}

sub command {
	my ( $self, $command ) = @_;
	if ( $self->{'DEBUG'} ) {
		print "$command\n";
	}
	else {
		system("echo $command");
		system("$command");
	}
}

sub get_time {
	my $dir = shift;
	opendir( DIR, $dir ) or print $!, "\n";
	my @matches = grep( /\d{10}/, readdir(DIR) );
	closedir(DIR);
	if(scalar @matches > 0){
		return $1 if $matches[0] =~ m/(\d{10})/;
	}
	return time;
}

sub get_all_ids{
	my ( $self ) = @_;
	my $dir = $self->ids_dir();
	my @files = <$dir/*>;
	my @ids;
	for my $file(@files){
		push (@ids, $self->get_sge_id($file));		
	}
	return \@ids;
}

sub read_intervals {
	my ( $self ) = @_;
	my $file = $self->{'CONFIG'}->{'GATKGENOMEBED'};
	my @data;
	open IN, $file or die "Can't open $file\n";
	while (<IN>) {
		chomp;
		next if m/GL/;
		push( @data, "$1" ) if m/(.+)\t(\d+)\t(\d+)/;
	}
	close IN;
	return @data;
}

return 1;
