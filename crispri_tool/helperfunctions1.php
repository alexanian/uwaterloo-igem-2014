<?php
//---------------------------------------------------
// This function finds the complementary strand to the sequence
function complement($input_seq) {
$step1 = str_replace("C","g",$input_seq);
$step2 = str_replace("G","c",$step1);
$step3 = str_replace("T","a",$step2);
$step4 = str_replace("A","t",$step3);
$output_seq = strtoupper($step4);
return strrev($output_seq);
}

//---------------------------------------------------
// This function determines if results will be found or not
function willfind($input_seq, $pam_seq, $sg_len) {

$comp_seq = complement($input_seq); // Find the coding strand

$pos1 = strpos($input_seq,$pam_seq); // Find the location of the first PAM sites, if they exist
$pos2 = strpos($comp_seq,$pam_seq);

if ($pos1 === false &&  $pos2 === false) {return false;} // If PAM does not exist in the sequence, do not proceed
else {return true;} // If it does, proceed
}

//---------------------------------------------------
// This function determines the form of the output of sgRNA...for now it's magic
function find_sg_seq($input_seq, $pos, $pam_seq, $sg_len) {
$startpos = $pos-$sg_len;
/*
if ($startpos < 0) {
	$end = $sg_len-$pos;
	$output_string = substr($input_seq, 0, $pos+2);
}
*/
if ($startpos < 0) {
	$output_string = "Not long enough to target :(";
}
else {
	$len = $sg_len;
	$output_string = substr($input_seq, $startpos+1, $len+strlen($pam_seq)-1);
}
return $output_string;
}


//---------------------------------------------------
// This function finds the locations of all PAM sites that may be targeted
function findsites($input_seq, $pam_seq, $sg_len) {

$comp_seq = complement($input_seq); // Find the coding strand

$i = 0; // Dummy index
$breaker = true; // Condition for when to exit loop
// Loop through the string
while ($breaker === true) {
	$pos1 = strpos($input_seq,$pam_seq,$i); // Finds position of PAM substring on template strand
	$pos2 = strpos($comp_seq,$pam_seq,$i); // Finds position of PAM substring on coding strand
	
	//--------------------------
	// Case 1: The PAM site is not found on either strand
	if ($pos1 === false && $pos2 === false) {
		$breaker = false;
	}
	
	//--------------------------
	// Case 2: PAM found only on coding strand
	elseif ($pos1 === false) {
		$i = $pos2+1;
		$sgRNA = find_sg_seq($comp_seq, $pos2, $pam_seq, $sg_len);
		echo "<tr> 
			<td>$i</td>
			<td>Coding</td>
			<td>$sgRNA</td>
			</tr>";
	}
	
	//--------------------------
	// Case 3: PAM found only on template strand
	elseif ($pos2 === false) {
		$i = $pos1+1;
		$sgRNA = find_sg_seq($input_seq, $pos1, $pam_seq, $sg_len);
		echo "<tr>
			<td>$i</td>
			<td>Template</td>
			<td>$sgRNA</td>
			</tr>";
	}
	
	//--------------------------
	// Case 4: Both still have PAM sites to find...
	//...but both occur at same location
	elseif ($pos1 === $pos2) {
		$i = $pos1+1;
		$sgRNA1 = find_sg_seq($input_seq, $pos1, $pam_seq, $sg_len);
		$sgRNA2 = find_sg_seq($comp_seq, $pos2, $pam_seq, $sg_len);
		echo "<tr> 
			<td>$i</td>
			<td>Template</td>
			<td>$sgRNA1</td>
			</tr>
			<tr> 
			<td>$i</td>
			<td>Coding</td>
			<td>$sgRNA2</td>
			</tr>";
	}
	//...else find which has a PAM occuring next and print that PAM
	else {
		// The next PAM is on the template strand
		if ($pos1 < $pos2) {
			$i = $pos1+1;
			$strand = "Template";
			$sgRNA = find_sg_seq($input_seq, $pos1, $pam_seq, $sg_len);
		}
		// The next PAM is on the coding strand
		else {
			$i = $pos2+1;
			$strand = "Coding";
			$sgRNA = find_sg_seq($comp_seq, $pos2, $pam_seq, $sg_len);
		}
		echo "<tr> 
			<td>$i</td>
			<td>$strand</td>
			<td>$sgRNA</td>
			</tr>";
				
	} //Matches else for case 4
} //Matches while loop

return true;

}

?>