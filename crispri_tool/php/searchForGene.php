<?php

require 'sql_helper.php';

function processRequest( $genename ) { 
    $sg_len = 20;
    
    $PossibleTargets = new PossibleTargetsDB();
    echo json_encode($PossibleTargets->searchGene($genename));
}

$genename = strtoupper($_POST['GeneName']);

processRequest( $genename );
?>
