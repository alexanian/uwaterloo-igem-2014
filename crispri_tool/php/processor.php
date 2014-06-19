<?php

require 'early_fxns.php';
require 'sql_helper.php';

function processRequest( $genenames, $genes, $promoters, $PAM ) { 
    $sg_len = 20;
    
    $PossibleTargets = new PossibleTargetsDB();
    $targetCache = new TargetCache();
	$a = 0;
	
    foreach ($genenames as $a => $b) {
	    $coding = strtoupper($promoters[$a] . $genes[$a]);
	    $num = findSitesUsingCodingStrand($PossibleTargets, $targetCache, $genenames[$a], $coding, $PAM, $sg_len);
    }

    echo json_encode($targetCache);
}

$PAM = strtoupper($_POST['pamseq']);
$genenames = $_POST['GENnames'];
$promoters = $_POST['PROs'];
$genes = $_POST['GENs'];

processRequest( $genenames, $genes, $promoters, $PAM );
?>
