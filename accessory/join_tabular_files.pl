use strict;
use Getopt::Long;
use GeneAnnotator;
use Data::Dumper;
####### get arguments      ###
my ( $table, $annotation, $table_id_columns, $table_columns,
	$annotation_id_columns, $annotation_columns, $skip_table_header, $skip_annotation_header, $annotation_header );
GetOptions(
	'table=s'                 => \$table,
	'annotation=s'            => \$annotation,
	'table_id_columns=s'      => \$table_id_columns,
	'annotation_id_columns=s' => \$annotation_id_columns,
	'annotation_columns=s'      => \$annotation_columns,
	'table_columns=s'           => \$table_columns,
	'skip_table_header'             => \$skip_table_header,
	'skip_annotation_header'             => \$skip_annotation_header,
	'annotation_header=s' => \$annotation_header,
);

my $info = {};
my @head_print;

open ANN, "<$annotation" or die "Can't open $annotation\n";
if ( !$skip_annotation_header ) {
	my $header = <ANN>;
	my $head_elements = get_elements( $header, $annotation_columns );
	@head_print = @$head_elements;
}
if($annotation_header){
	@head_print = split(/,/, $annotation_header);
}
while (<ANN>) {
	my $id  = get_id( $_, $annotation_id_columns );
	my $ann = get_elements(  $_, $annotation_columns );
	$info->{ $id } = $ann;
}
close ANN;
open TABLE, "<$table" or die "Can't open $table\n";
if ( !$skip_table_header ) {
	my $header = <TABLE>;
	my $head_elements = get_elements( $header, $table_columns );
	@head_print = (@$head_elements,@head_print) if $head_elements;
	print join( "\t", @head_print ), if scalar @head_print;
	print "\n";
}
while (<TABLE>) {
	my $id  = get_id( $_, $table_id_columns );
	my $table = get_elements(  $_, $table_columns );
	my @to_print = (@$table);
	my @to_print = (@$table, @{$info->{$id}}) if $info->{$id};
	print join( "\t", (@to_print)), "\n";
}

close TABLE;

sub get_id {
	my ( $array_str, $columns ) = @_;
	return join( '_', @{get_elements( $array_str, $columns )} );
}

sub get_elements {
	my ( $array_str, $columns )  = @_;
	chomp $array_str;
	my @array = split /\t/, $array_str;
	my @nums  = split /,/,  $columns;
	my @data = map { $array[$_] } @nums;
	return \@data;
}
