use strict;
use Data::Dumper;
use XML::DOM;

my $file       = "uniprot_sprot.xml";
my $dom_parser = new XML::DOM::Parser;
my $doc        = $dom_parser->parsefile($file);
my $root       = $doc->getDocumentElement();
my $entries    = $doc->getElementsByTagName("entry");

my @list = (
	'accession',          'gene',
	'fullName',           'function',
	'catalytic activity', 'tissue specificity',
	'pathway',            'disease'
);
print join( "\t", @list ), "\n";
for ( my $e = 0 ; $e < $entries->getLength ; $e++ ) {
	my $entry = $entries->item($e);

	my $org =
	  $entry->getElementsByTagName("organism")->item(0)
	  ->getElementsByTagName("name")->item(0)->getFirstChild()->getData;
	next unless $org eq 'Homo sapiens';  
	#accession
	my $accession =
	  $entry->getElementsByTagName("accession")->item(0)->getFirstChild()
	  ->getData;

	#fullName
	#uniprot->entry->protein->recommendedName->fullName
	my $fullName =
	  $entry->getElementsByTagName("fullName")->item(0)->getFirstChild()
	  ->getData;

	#gene
	#uniprot->entry->gene->name
	my $gene =
	  $entry->getElementsByTagName("gene")->item(0)
	  ->getElementsByTagName("name")->item(0)->getFirstChild()->getData;

	#All comments
	my $comments = $entry->getElementsByTagName("comment");
	my $types    = {
		'function'           => [],
		'catalytic activity' => [],
		'tissue specificity' => [],
		'pathway'            => [],
		'disease'            => [],
		'gene'               => [],
		'fullName'           => [],
		'accession'          => [],
		''                   => [],
	};
	push( @{ $types->{'accession'} }, $accession );
	push( @{ $types->{'gene'} },      $gene );
	push( @{ $types->{'fullName'} },  $fullName );
	for ( my $i = 0 ; $i < $comments->getLength ; $i++ ) {
		my $comment = $comments->item($i);
		my $type    = $comment->getAttribute('type');
		my $value   = $comment->getElementsByTagName("text")->item(0);
		if ($value) {
			my $value_data = $value->getFirstChild()->getData;
			push( @{ $types->{$type} }, $value_data );
		}
	}
	to_string($types);

	#function
	#uniprot->entry->comment type="function"->text
	#'catalytic activity'
	#'tissue specificity'
	#'pathway',
	#'disease'
}

sub to_string {
	my ($data) = @_;
	my @to_print;
	for my $el (@list) {
		push( @to_print, join( ",", @{ $data->{$el} } ) );
	}
	print join( "\t", @to_print ), "\n";
}

