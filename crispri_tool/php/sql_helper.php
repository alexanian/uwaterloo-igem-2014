<?php

class PossibleTargetsDB
{
    private $SQLUser    = 'root';
    private $SQLDb      = 'possible_targets';
    private $SQLPass    = '';
    private $SQLHost    = 'localhost';
    private $SQLCon     = null;
    
    /**
     * Construct PossibleTargets DB Helper.
     * Will kill script ('die') if unable to connect.
     */
    function __construct() {
        $this->SQLCon = new mysqli( $this->SQLHost, $this->SQLUser, $this->SQLPass, $this->SQLDb );
        
        if($this->SQLCon->connect_errno) {
            die('Failed to connect to SQL database');   
        }
    }
    
    public function searchGene( $gene ) { 
        $geneEscaped = str_replace(['%', '_'], ['\%','\_'], $gene);
        $findGeneQuery = "SELECT GeneId, GeneName, CodingSeq AS CodingSequence FROM GeneTable WHERE GeneName LIKE '%$geneEscaped%'";
        $geneList = [];
        if( $result = $this->SQLCon->query($findGeneQuery) ) { 
            $geneList = $result->fetch_all(MYSQLI_ASSOC);
            $result->close();
            
            $i = 0;
            foreach($geneList as $i => $b) { 
                $geneList[$i]['GeneId'] = intval($b['GeneId']);
            }
        }
        return $geneList;
    }
    
    public function queryTargets( $geneid, $sg_len ) { 
        $targetQuery = "SELECT GeneTable.GeneId, GeneName, Strand, Position, Result AS sgRNA, Repression
                        FROM Targets
                        JOIN GeneTable ON GeneTable.GeneId = Targets.GeneId
                        JOIN StrandTable ON StrandTable.StrandId = Targets.StrandId
                        WHERE   Targets.GeneId = $geneid
                            AND (Position - $sg_len) >= 0";
        $targetList = [];
        if( $result = $this->SQLCon->query($targetQuery) ) { 
            $targetList = $result->fetch_all(MYSQLI_ASSOC);
            $result->close();
            
            $i = 0;
            foreach($targetList as $i => $b) { 
                $targetList[$i]['GeneId'] = intval($b['GeneId']);
                $targetList[$i]['Position'] = intval($b['Position']);
                $targetList[$i]['Repression'] = doubleval($b['Repression']);
            }
        }
        return $targetList;
    }
    
    public function insertGene( $gene, $gene_code_seq, &$geneId ) {
        $geneFound = null;
        $geneId = 0;
        // Find out if it's already in the database
        foreach( $this->searchGene( $gene ) as $possibleGene ) { 
            if( $possibleGene['CodingSequence'] == $gene_code_seq &&
                $possibleGene['GeneName'] == $gene ) {
                // Just to be sure we have found the right gene. We really shouldn't need this.
                // @TODO : Front-end does a lookup first and determines conflicts, then picks an id or doesn't;
                $geneFound = $possibleGene;
                $geneId = $geneFound['GeneId'];
                return true;
            }
        }
        
        if($geneId == 0) {
            // Aww...we didn't find it...time to do work.
            $insertGeneQuery = "INSERT INTO GeneTable ( GeneName, CodingSeq ) VALUES ( '$gene', '$gene_code_seq' ); SELECT LAST_INSERT_ID() AS GeneId;";
            if( $this->SQLCon->multi_query($insertGeneQuery) ) {
                $this->SQLCon->next_result();
                $result = $this->SQLCon->store_result();
                $obj = $result->fetch_object();
                $result->close();
                
                $geneId = intval($obj->GeneId);
            } else { 
                die( "Unable to insert gene '$gene' into database." );
            }
        }
        return false;
    }
    
    /**
     * Inserts a Gene Target Record into the database
     * @param   geneid  The Gene Id in Database
     * @param   strand  The strand the target is on (0 if Coding, 1 if template)
     * @param   pos     Position to target
     * @param   sgRNA   The targeting sgRNA sequence as string
     * @returns bool    True if successful, False otherwise.
     */
    public function insertGeneTarget($geneid, $strand, $pos, $sgRNA) {
        $strandId = $strand + 1;
        $insertRecord = "INSERT INTO Targets (GeneId, StrandId, Position, Result)
                        VALUES ('$geneid', '$strandId', '$pos', '$sgRNA')";
        if( !$this->SQLCon->real_query($insertRecord) ) { 
            die( "Failed to insert target into database for gene id '$geneid'" );
        }
    }
}

?>
