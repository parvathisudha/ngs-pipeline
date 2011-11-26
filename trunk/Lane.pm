package Lane;
use strict;

sub new {
	my ( $class, $dir, $num, $archived, $paired, $rg ) = @_;

	my $self = {
		DIR            => $dir,
		NUM            => $num,
		ARCHIVED       => $archived,
		PAIRED       => $paired,
		RG       => $rg,
	};
	bless $self, $class;
	return $self;
}

return 1;
