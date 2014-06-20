<?php

class TargetCache
{
    public  $data               = [];
    public  $recordsTotal       = 0;
    public  $recordsFiltered    = 0;
    public  $draw               = 1;
    
    public function addTargets( $targetList ) { 
        foreach( $targetList as $target ) { 
            $this->addTargetObject( $target );
        }
    }
    
    public function addTargetObject( $target ) { 
        $this->recordsTotal++;
        $this->recordsFiltered++;
        
        $this->data[] = $target;
    }
    
    public function addTarget( $gene, $strand, $position, $sgRNA ) { 
        $this->recordsTotal++;
        $this->recordsFiltered++;
        
        $this->data[] =
            array(
                "GeneName" => $gene, 
                "Strand" => $strand,
                "Position" => $position,
                "sgRNA" => $sgRNA,
                );
    }
}

//---------------------------------------------------
// This function finds the complementary strand to the sequence
function complement($input_seq) {
    $step1 = str_replace("C","G",$input_seq);
    $step2 = str_replace("G","C",$step1);
    $step3 = str_replace("T","A",$step2);
    $step4 = str_replace("A","T",$step3);
    $output_seq = strtoupper($step4);
    return strrev($output_seq);
}

/**
 * This function determines the form of the output of sgRNA...for now it's magic
 * @param   input_seq           Input Sequence
 * @param   pos                 Position of Site
 * @param   pam_seq             PAM Sequence
 * @param   sg_len              Minimum length of sgRNA
 * @returns array               { sgRNA => sgRNA Sequence, valid => If it is of valid length }
 */
function find_sg_seq($input_seq, $pos, $pam_seq, $sg_len) {

    $startpos = $pos - $sg_len;
    $valid = ($startpos >= 0);
    
    if ($valid) {
	    $len = $sg_len;
	    $output_string = substr($input_seq, $startpos+1, $len+strlen($pam_seq)-1);
    }

    return array( 
        "sgRNA" => $output_string,
        "valid" => $valid
        );
}

/**
 * This function finds the locations of all PAM sites that may be targeted
 * @param   possibleTargets     PossibleTargetsDB instance
 * @param   targetCache         Target Cache instance
 * @param   geneid              Gene Id in Database
 * @param   genename            Gene Name
 * @param   original_seq        Coding Strand Sequence
 * @param   strand              0 if Coding, 1 if Template
 * @param   pam_seq             PAM Sequence
 * @param   sg_len              Minimum sgRNA Length
 */
function findsites(&$possibleTargets, &$targetCache, $geneid, $genename, $original_seq, $strand, $pam_seq, $sg_len) {
    $i = 0; // Dummy index
    $num_found = 0;
    $pos = true;
    $input_seq = ($strand == 1 ? complement($original_seq) : $original_seq);
    $strandName = ($strand == 1 ? "Template" : "Coding");
    
    // Loop through the string
    // and Find position of PAM substring on strand
    // Break if nothing is found
	while( ($pos = strpos($input_seq, $pam_seq, $i)) !== false ) {
        
        // There is a target
        $i = $pos + 1;
	    $sgRNA = find_sg_seq($input_seq, $pos, $pam_seq, $sg_len);
	    
	    if( $sgRNA['valid'] ) {
	        $targetCache->addTargetObject( array( "GeneName" => $genename, "Strand" => $strandName, "Position" => $pos, "sgRNA" => $sgRNA['sgRNA'] ) );
	    }
	
	    $possibleTargets->insertGeneTarget($geneid,$strand,$i,$sgRNA['sgRNA']);
    }
}

/**
 * This function finds the locations of all PAM sites that may be targeted on a Gene
 * It will first lookup the gene in the database and see if the data is already cached.
 * If it is found, it will cache any valid targets for returning as json in the TargetCache.
 * Otherwise, it will find sites, insert them into the DB, and cache the valid ones in the targetcache.
 * @param   possibleTargets     PossibleTargetsDB instance
 * @param   targetCache         Target Cache instance
 * @param   genename            Gene Name
 * @param   original_seq        Coding Strand Sequence
 * @param   pam_seq             PAM Sequence
 * @param   sg_len              Minimum sgRNA Length
 */
function findSitesUsingCodingStrand(&$possibleTargets, &$targetCache, $genename, $coding, $pam_seq, $sg_len) {
    $geneId = 0;
    
    if( $possibleTargets->insertGene( $genename, $coding, $geneId ) ) { 
        $targetCache->addTargets( $possibleTargets->queryTargets( $geneId, $sg_len ) );
    } else { 
        findsites($possibleTargets, $targetCache, $geneId, $genename, $coding, 0, $pam_seq, $sg_len);
        findsites($possibleTargets, $targetCache, $geneId, $genename, $coding, 1, $pam_seq, $sg_len);
    }
}
?>
