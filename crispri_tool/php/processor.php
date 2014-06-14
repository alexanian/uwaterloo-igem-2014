<?php

require 'early_fxns.php';
require 'sql_helper.php';

$PAM = strtoupper($_POST['pamseq']);
$sg_len = 20;

$genenames = $_POST['GENnames'];
$promoters = $_POST['PROs'];
$genes = $_POST['GENs'];
$a=0;

/*
echo "The PAM site is " . $PAM . "<br>";
echo "The gene names are " . $genenames[$a] . "<br>";
echo "The promoters are " . $promoters[$a] . "<br>";
echo "The genes " . $genes[$a] . "<br>";
*/

connect_db();

echo "<table id='resultsTable' class='display'>
	<thead>
		<tr>
			<th>Gene</th>
			<th>Strand</th>
			<th>Position</th>
			<th>sgRNA</th>
			<th>Effectiveness</th>
		</tr>
	</thead>
	<tbody>";
	
foreach ($genenames as $a => $b) {
	$coding = strtoupper($promoters[$a] . $genes[$a]);
	$template = complement($coding);
	$num1 = findsites($genenames[$a],$coding,"Coding",$PAM,$sg_len);
	echo "The coding strand of $genenames[$a] has been processed with $num1 possible targets found.<br>";
	$num2 = findsites($genenames[$a],$template,"Template",$PAM,$sg_len);
	echo "The template strand of $genenames[$a] has been processed with $num2 possible targets found.<br><br>";	
}

echo "</tbody>
	</table>
	
	<script type='text/javascript' language='javascript' class='init'>
	$(document).ready(function() {
		$('#resultsTable').dataTable();
	} );
	</script>";

?>