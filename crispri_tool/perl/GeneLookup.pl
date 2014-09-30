use strict;
use warnings;

use Bio::EnsEMBL::LookUp;

my $lookup = Bio::EnsEMBL::LookUp->new();
my $dbas = $lookup->get_by_name_exact(@ARGV[0]);
my @genes = @{$dbas->get_GeneAdaptor()->fetch_all_by_description('%' . @ARGV[1] . '%')};

print "[\n";
my $c = 0;
Genes: while( scalar(@genes) ) { 
    print "{ id:\"" . @genes[$c]->stable_id . "\", description:\"" . @genes[$c]->description . "\", sequence:\"" . @genes[$c]->seq . "\"}";
    $c++;
  if($c >= scalar(@genes) ) { last Genes; }
    print ",\n"
}
print "\n]";
