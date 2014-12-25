package Omim;
use strict;
use Data::Dumper;
use Record;
use Field;
use Morbidmap;

sub new {
	my ( $class, $dir) = @_;
	my $self = {
		dir => $dir,
		omim2data => {},
		gene2omim => {},
	};
	bless $self, $class;
	$self->read_morbidmap();
	$self->read_omim();
	return $self;
}

#14.325|2|14|11|14q24.2|PSEN1, AD3|C|Presenilin 1||104311|Fd|||Alzheimer disease, type 3, 607822 (3); Alzheimer disease, type 3,|with spastic paraparesis and unusual plaques, 607822 (3); Alzheimer disease, type 3, with spastic paraparesis and apraxia, 607822 (3); Dementia,|frontotemporal, 600274 (3); Pick disease, 172700 (3); Cardiomyopathy, dilated, 1U, 613694 (3); Acne inversa, familial, 3, 613737 (3)||

sub get_morbid_record_by_gene_name{
	my ($self, $gene_name) = @_;
	my $record_id = $self->{gene2omim}->{$gene_name};
	if($record_id){
		return $self->{omim2data}->{$record_id};
	}
	else{
		return undef;
	}
}

sub read_morbidmap{
        my ($self) = @_;
        my $dir = $self->{dir};
        my $morbidmap_file = "$dir/morbidmap";
        open MORBIDMAP, "<$morbidmap_file" or die "Can't open omim.txt in $morbidmap_file\n";
        while(my $str = <MORBIDMAP>){
		my $morbid = Morbidmap->new($str);
		$self->{gene2omim}->{$morbid->{gene_name}} = $morbid->{phenotype_id};
	}
	close MORBIDMAP;
}

sub read_omim{
	my ( $self) = @_;
	my $dir = $self->{dir};
	my $omimtxt_file = "$dir/omim.txt";
	open OMIM, "<$omimtxt_file" or die "Can't open omim.txt in $omimtxt_file\n";
	my $record_buffer = "";
	while(my $str = <OMIM>){
		if($str =~ m/\*RECORD\*/){
			if($record_buffer){
				my $omim_record = Record->new($record_buffer);
				my $id = $omim_record->get_id();
				$self->{omim2data}->{$id} = $omim_record;
				#chomp $str;
				$record_buffer = "";
			}
		}
		else{
			$record_buffer .= $str;
		}
	}
	my $omim_record = Record->new($record_buffer);
	my $id = $omim_record->get_id();
	$self->{omim2data}->{$id} = $omim_record;
	close OMIM;
}


1;
