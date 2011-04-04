package Project;
use strict;
#file names for all tasks

sub new {
	my ( $class, $config, $debug) = @_;
	my $time = time;
	my $self = {
		CONFIG => $config,
		DEBUG => $debug,
		TIME => $time,
	};
	bless $self, $class;
	return $self;
}

#folders

sub dir{
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'};
}
sub script_dir{
	my ($self) = @_;
	return $self->dir() . '/tasks';
}
sub output_dir{
	my ($self) = @_;
	return $self->dir() . '/out';
}
sub error_dir{
	my ($self) = @_;
	return $self->dir() . '/errors';
}
sub ids_dir{
	my ($self) = @_;
	return $self->dir() . '/ids';
}

sub task_id_file{
	my ( $self, $job_name ) = @_;
	return $self->ids_dir() . "/$job_name.id";
}

sub task_id{
	my ( $self, $job_name ) = @_;
	my $file = $self->task_id_file($job_name);
	return $self->get_sge_id($file);
}

sub make_folder_structure{
	my ($self) = @_;
	$self->mkdir($self->dir());
	$self->mkdir($self->script_dir());
	$self->mkdir($self->output_dir());
	$self->mkdir($self->error_dir());
	$self->mkdir($self->ids_dir());
}

sub mkdir{
	my ($self, $directory) = @_;
	$self->command("mkdir $directory");
}

sub forward_name{
	my ($self, $lane) = @_;
	return $lane->forward_name();
}

sub reverse_name{
	my ($self, $lane) = @_;
	return $lane->reverse_name();
}

sub forward_prefix{
	my ($self, $lane) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->forward_align_id($lane);
}

sub reverse_prefix{
	my ($self, $lane) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->reverse_align_id($lane);
}

sub forward_sai{
	my ($self, $lane) = @_;
	return $self->forward_prefix($lane) . ".sai";
}

sub reverse_sai{
	my ($self, $lane) = @_;
	return $self->reverse_prefix($lane) . ".sai";
}

sub get_read_id{
	my ( $self, $lane, $is_forward ) = @_;
	my $read_direction = $is_forward ? 1 : 2;
	#s_$num_$read_direction_sequence.txt
	my $name =
	    $self->_get_id($lane->{DIR}) . '_' . 's_'
	  . $lane->{NUM} . '_'
	  . $read_direction
	  . '_sequence.txt';
	return $name;
}

sub forward_align_id{
	my ($self, $lane) = @_;
	my $id = 'aln.' . $self->_get_id($self->get_read_id($lane, 1));
	return $id;
}

sub reverse_align_id{
	my ($self, $lane) = @_;
	return 'aln.' . $self->_get_id($self->get_read_id($lane, 0));
}

sub sam{
	my ($self, $lane) = @_;
	return $self->forward_prefix($lane) . ".sam";
}

sub bam{
	my ($self, $lane) = @_;
	return $self->forward_prefix($lane) . ".bam";
}

sub sorted_prefix{
	my ($self, $lane) = @_;
	return $self->forward_prefix($lane) . ".sorted";
}

sub sorted{
	my ($self, $lane) = @_;
	return $self->sorted_prefix($lane) . ".bam";
}

sub file_prefix{
	my ($self) = @_;
	return $self->{'CONFIG'}->{'DIR'} . '/' . $self->{'CONFIG'}->{'PROJECT'};
}

sub merged{
	my ($self) = @_;
	return $self->file_prefix() . ".bam";
}

sub merged_id{
	my ($self) = @_;
	return 'merge.' . $self->_get_id($self->merged());
}
sub merged_indexed{
	my ($self) = @_;
	return $self->file_prefix() . ".bai";
}
sub merged_indexed_id{
	my ($self) = @_;
	return 'index_merged.' . $self->_get_id($self->merged_indexed());
}

sub merged_sorted_prefix{
	my ($self) = @_;
	return $self->file_prefix() . ".sorted";
}

sub merged_sorted{
	my ($self) = @_;
	return $self->merged_sorted_prefix() . ".bam";
}
sub merged_sorted_id{
	my ($self) = @_;
	return 'sort_merged.' . $self->_get_id($self->merged_sorted());
}

sub gatk_vcf{
	my ($self) = @_;
	return $self->file_prefix() . ".vcf";
}

sub gatk_vcf_id{
	my ($self) = @_;
	return 'gatk_snps.' . $self->_get_id($self->gatk_vcf());
}

sub depth_coverage{
	my ($self) = @_;
	return $self->file_prefix() . ".depth_coverage";
}

sub depth_coverage_id{
	my ($self) = @_;
	return 'depth_cov.' . $self->_get_id($self->depth_coverage());
}


sub callable_loci{
	my ($self) = @_;
	return $self->file_prefix() . ".callable_loci";
}
sub callable_loci_summary{
	my ($self) = @_;
	return $self->file_prefix() . ".callable_loci_summary";
}

sub callable_loci_id{
	my ($self) = @_;
	return 'callable.' . $self->_get_id($self->callable_loci());
}

sub eff_vcf{
	my ($self) = @_;
	return $self->file_prefix() . ".eff.vcf";
}

sub eff_vcf_id{
	my ($self) = @_;
	return 'effect.' . $self->_get_id($self->eff_vcf());
}

sub filter_snps{
	my ($self) = @_;
	return $self->file_prefix() . ".eff.filtered.vcf";
}

sub filter_snps_id{
	my ($self) = @_;
	return 'filter.' . $self->_get_id($self->filter_snps());
}

sub genome_coverage{
	my ($self) = @_;
	return $self->file_prefix() . ".genome_coverage.bed";
}

sub genome_coverage_id{
	my ($self) = @_;
	return 'g_cov.' . $self->_get_id($self->genome_coverage());
}

sub genome_coverage_bga{
	my ($self) = @_;
	return $self->file_prefix() . ".genome_coverage_bga.bed";
}

sub genome_coverage_bga_id{
	my ($self) = @_;
	return 'g_cov_bga.' . $self->_get_id($self->genome_coverage_bga());
}

sub move_bedtools_results_id{
	my ($self) = @_;
	return 'mv_bed_res.' . $self->_get_id('move_bedtools_results_id');	
}

sub all_indexed_ids{
	my ($self) = @_;
	my $lanes = $self->{'CONFIG'}->{'LANES'};
	my @ids = "";
	for my $lane(@$lanes){
		push (@ids, $self->task_id($self->index_id($lane)));
	}
	return join (',', @ids);
}

sub sorted_id{
	my ($self, $lane) = @_;
	return 'sort.' . $self->_get_id($self->sorted($lane));
}

sub sam_id{
	my ($self, $lane) = @_;
	return 'sam.' . $self->_get_id($self->sam($lane));
}

sub import_id{
	my ($self, $lane) = @_;
	return 'import.' . $self->_get_id($self->bam($lane));
}
sub bai{
	my ($self, $lane) = @_;
	return $lane->forward_name() . ".sorted.bam.bai";
}
sub index_id{
	my ($self, $lane) = @_;
	return 'index.' . $self->_get_id($self->bai($lane));
}

sub get_garbage_files{
	my ($self) = @_;
	my @array;
	my $lanes = $self->{'CONFIG'}->{'LANES'};
	for my $lane(@$lanes){
		push(@array, $self->forward_name($lane));
		push(@array, $self->reverse_name($lane));
		push(@array, $self->sam($lane));
		push(@array, $self->bam($lane));
		push(@array, $self->sorted($lane));
		push(@array, $self->bai($lane));
		
	}
	push(@array, $self->merged());		
	return \@array;
}

sub clean_id{
	my ($self) = @_;
	return 'clean.' . $self->_get_id('clean_id');	
}



sub _get_id{
	my ($self, $file) = @_;
	$file =~ s/\//_/g;
	$file =~ s/\://g;
	$file .= "_" . $self->{'TIME'};
	return $file;
}
sub get_sge_id{
	my ($self, $file) = @_;
	if($self->{'DEBUG'}){
		$file = "id.txt";
	}
	open IN, "<$file" or die "Can't open file: $file for reading task id!\n";
	while(<IN>){
		return $1 if m/Your job\s(\d+)\s/;
	}
	close IN;
}

sub command{
        my ($self, $command) = @_;
        if($self->{'DEBUG'}){
        	print "$command\n";
        }        
        else{
	        system("echo $command");
	        system("$command");       	
        }
}

return 1;