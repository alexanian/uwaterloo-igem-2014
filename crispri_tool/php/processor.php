<?php

require 'find_targets_fxns.php';

function processSeqSearch( $genenames, $genes, $promoters, $PAM ) { 
    $sg_len = 20;
	
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

processSeqSearch( $genenames, $genes, $promoters, $PAM );
?>
