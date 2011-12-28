use XML::Simple;
use Data::Dumper;
my $config = XMLin("S000043.20111130.xml");

print Dumper($config);
