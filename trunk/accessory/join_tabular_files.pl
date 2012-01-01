use strict;
use Getopt::Long;
use GeneAnnotator;
use Data::Dumper;
####### get arguments      ###
my ( $table, $annotation, $table_id_columns, $table_columns,
	$annotation_id_columns, $annotation_columns, $skip_header );
GetOptions(
	'table=s'                 => \$table,
	'annotation=s'            => \$annotation,
	'table_id_columns=s'      => \$table_id_columns,
	'annotation_id_columns=s' => \$annotation_id_columns,
	'annotation_columns=s'      => \$annotation_columns,
	'table_columns=s'           => \$table_columns,
	'skip_header'             => \$skip_header,
);

my $info = {};
my @head_print;

open ANN, "<$annotation" or die "Can't open $annotation\n";
if ( !$skip_header ) {
	my $header = <ANN>;
	my $head_elements = get_elements( $header, $annotation_columns );
	@head_print = @$head_elements;
}
while (<ANN>) {
	my $id  = get_id( $_, $annotation_id_columns );
	my $ann = get_elements(  $_, $annotation_columns );
	$info->{ $id } = $ann;
}
close ANN;
open TABLE, "<$table" or die "Can't open $table\n";
if ( !$skip_header ) {
	my $header = <TABLE>;
	my $head_elements = get_elements( $header, $table_columns );
	@head_print = (@head_print,@$head_elements) if $head_elements;
	print join( "\t", @head_print ), if scalar @head_print;
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
