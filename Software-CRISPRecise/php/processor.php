<?php
require 'find_targets_fxns.php';

function processSeqSearch( $gene, $geneSeq, $promoterSeq, $PAM ) { 
    $sg_len = 20;
	
	echo "$gene";
	echo "$PAM";
	echo "$promoterSeq";
	
	$targetCache = new TargetCache();
	
	$num = findSitesUsingCodingStrand($targetCache, $gene, $geneSeq, $promoterSeq, $PAM, $sg_len);
}

$PAM = strtoupper($_POST['pamseq']);
$gene = $_POST['Gene'];
$promoterSeq = $_POST['PromoterSeq'];
$geneSeq = $_POST['GeneSeq'];

processSeqSearch($gene, $geneSeq, $promoterSeq, $PAM );
?>
