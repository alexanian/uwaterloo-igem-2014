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

/**
 * This function finds the reverse complement strand to the sequence
 * @param	input_seq			Input Sequence
 * @returns	output_seq			Reverse Complement of input_seq
 */
function complement($input_seq) {
    $complements = array("C" => "G", "G" => "C", "T" => "A", "A" => "T");
    $output_seq = strtoupper(strtr($input_seq, $complements));
    return strrev($output_seq);
}
/**
 * This function returns a pattern corresponding to possibly ambiguous
 * nucleotide input, e.g. NGG -> .GG
 * @param	PAM					PAM Sequence (may contain ambiguous characters)
 * @returns	pam_seq				Pattern that will match all appropriate PAM sequences
 */
function nucleotide2pattern($PAM) {
	$patterns = array("W" => "[AT]", "S" => "[CG]",
					"M" => "[AC]", "K" => "[GT]",
					"R" => "[AG]", "Y" => "[CT"],
					"B" => "[CGT]", "D" => "[AGT]",
					"H" => "[ACT]", "V" => "[ACG]",
					"N" => ".", "-" => ".");
	$pam_seq = strtr($PAM, $patterns);
	return($pam_seq);	
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
 * @param   targetCache         Target Cache instance
 * @param   geneid              Gene Id in Database
 * @param   genename            Gene Name
 * @param   original_seq        Coding Strand Sequence
 * @param   strand              0 if Coding, 1 if Template
 * @param   pam_seq             PAM Sequence
 * @param   sg_len              Minimum sgRNA Length
 */
function findsites(&$targetCache, $gene, $geneSeq, $promoterSeq, $strand, $PAM, $sg_len) {

	$strandName = ($strand == 1 ? "Template" : "Coding");
	
	// WARN: Currently not searching promoter sequence!
    $input_seq = ($strand == 1 ? complement($geneSeq) : $geneSeq); 

    // Convert ambiguous nucleotide characters (e.g. 'W') to regex
    $pam_seq = nucleotide2pattern($PAM);	
    
	// Locate all matches to PAM sequence
	$foundMatch = preg_match($pam_seq, $input_seq, $matches, PREG_OFFSET_CAPTURE);

	if($foundMatch) {
		//There is a target
		foreach ($matches[0] as $position) {
			$sgRNA = find_sg_seq($input_seq, $position, $pam_seq, $sg_len);
		}
		if( $sgRNA['valid'] ) {
	        $targetCache->addTargetObject( array( "GeneName" => $genename,
			"Strand" => $strandName, "Position" => $pos, "sgRNA" => $sgRNA['sgRNA'] ) );
	    }
	}
}
/**
 * This function finds the locations of all PAM sites that may be targeted on a Gene
 * It will first lookup the gene in the database and see if the data is already cached.
 * If it is found, it will cache any valid targets for returning as json in the TargetCache.
 * Otherwise, it will find sites, insert them into the DB, and cache the valid ones in TargetCache.
 * @param   targetCache         Target Cache instance
 * @param   gene            Gene Name
 * @param   geneSeq	        Coding Strand Gene Sequence
 * @param   promoterSeq		Coding Strand Promoter Sequence
 * @param   PAM         	PAM Sequence
 * @param   sg_len          Minimum sgRNA Length
 */
function findSitesUsingCodingStrand(&$targetCache, $gene, $geneSeq, $promoterSeq, $PAM, $sg_len) {
	// The first search will find matches in the coding strand, while the second search will
    // take the reverse complement 
    findsites($gene, $geneSeq, $promoterSeq, $coding, 0, $PAM, $sg_len);
    findsites($gene, $geneSeq, $promoterSeq, 1, $PAM, $sg_len);
}
?>
