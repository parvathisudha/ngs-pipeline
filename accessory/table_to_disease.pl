use strict;
use Omim;
use Data::Dumper;
use Spreadsheet::WriteExcel;                             # Step 0
use Getopt::Long;

####### get arguments      ###
my ( $in, $out, $omim_dir );
GetOptions(
	'in=s' => \$in,
	'out=s' => \$out,
	'omim_dir=s' => \$omim_dir,
);

###### read OMIM data      ###
my $omim = Omim->new( "$omim_dir");

###### setup Excel file    ###
my $excel_file = $out;
system("rm -f $excel_file");
my $workbook = Spreadsheet::WriteExcel->new($excel_file); # Step 1
my $worksheet = $workbook->add_worksheet();               # Step 2
my $line_counter = 1;

##### add annotation to input file ###
open IN, "<$in" or die "Can't open $in\n";

my @display = qw/TI,no_header TX,description CS,inheritance/;

##### write header to excel file ###
my $in_arr_str = <IN>;
chomp $in_arr_str;
my @in_arr = split ("\t", $in_arr_str);
my @header = (@in_arr,@display);
my $first_cell = get_cell();
$worksheet->write_row($first_cell, \@header);

my $gene_header_num = get_column_for_header(\@in_arr,"HGNC");

##### write lines with annotations
while(my $curr_str = <IN>){
	chomp $curr_str;
	my @curr_arr = split ("\t", $curr_str);
	my $gene = $curr_arr[$gene_header_num];
	my $morbid = $omim->get_morbid_record_by_gene_name($gene);
	my @data;
	for( my $i = 0; $i < scalar @display; $i++){
		if($morbid){
			push(@data,$morbid->get_field($display[$i]));
		}
		else{
			push(@data,"");
		}
	}
	my $first_cell = get_cell();
	my @line_data = (@curr_arr,@data);
        $worksheet->write_row($first_cell, \@line_data);
}

close IN;

################################################
#	Subroutines
################################################
sub get_column_for_header{
	my ($arr, $header) = @_;
	for (my $i = 0; $i < scalar @$arr; $i++){
		if($$arr[$i] eq $header){
			return $i;
		}
	}
	return undef;
}

sub get_cell{
	my $cell = "A" . $line_counter;
	$line_counter++;
	return $cell;
}
