use strict;
use warnings;

use Bio::EnsEMBL::LookUp;

my $lookup = Bio::EnsEMBL::LookUp->new();
my $temp = 'staphylococcus_aureus_.*';
my @dbas = @{$lookup->get_all_by_name_pattern(@ARGV[0] . ".*")};
my $limit = 0 + @ARGV[1];

print "[\n";
my $c = 0;
Species: while( scalar(@dbas) ) { 
    print "{ name:\"" . @dbas[$c]->species() . "\", id:\"" . @dbas[$c]->species_id . "\"}";
    $c++;
  if( ($c >= scalar(@dbas) ) or ( $c >= $limit ) ) { last Species; }
    print ",\n"
}
print "\n]";
