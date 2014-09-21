use strict;
use warnings;

use Bio::EnsEMBL::LookUp;

my $lookup = Bio::EnsEMBL::LookUp->new();
my $temp = 'staphylococcus_aureus_.*';
my @dbas = @{$lookup->get_all_by_name_pattern(@ARGV[0] . ".*")};
my $limit = 0 + @ARGV[1];

print "{\n";
for( my $c = 0; ( $c < scalar(@dbas) ) and ( $c < $limit ); $c++ ) { 
    print "{ name:\"" . @dbas[$c]->species() . "\", id:\"" . @dbas[$c]->species_id . "\"}\n";
}
print "}";
