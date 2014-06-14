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
// This function determines the form of the output of sgRNA...for now it's magic
function find_sg_seq($input_seq, $pos, $pam_seq, $sg_len) {

$startpos = $pos-$sg_len;

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
function findsites($genename, $input_seq, $strand, $pam_seq, $sg_len) {

$i = 0; // Dummy index
$num_found = 0;
$breaker = true; // Condition for when to exit loop
// Loop through the string
while ($breaker === true) {
	$pos = strpos($input_seq,$pam_seq,$i); // Finds position of PAM substring on strand
	
	// Case 1: No more targets
	if ($pos === false) {
		$breaker = false;
	}
	
	// Case 2: More targets
	else {
		$i = $pos+1;
		$sgRNA = find_sg_seq($input_seq, $pos, $pam_seq, $sg_len);
		if ($sgRNA !== "Not long enough to target :(") {$num_found++;}
		
		insert_row($genename,$strand,$i,$sgRNA);
		putout_row($genename,$strand,$i,$sgRNA);
	}
}

return $num_found;

}

function putout_row($genename,$strand,$i,$sgRNA) {
	
	echo "<tr>
			<td>$genename</td>
			<td>$strand</td>
			<td>$i</td>
			<td>$sgRNA</td>
			<td>0</td>
		</tr>";
	
}

?>