<html>

<head>
<title>Results</title>
</head>


<body>

<header>
	<h2>CRISPRi Design Tool: Results</h2>
</header>


<?php

require 'helperfunctions1.php';
// NOT YET - include 'sqlhelper.php';

$PAMseq = strtoupper($_POST["PAMin"]); // Retrieve the PAM sequence
$PROseq = strtoupper($_POST["PROin"]); // Retrieve the promoter sequence
$GENseq = strtoupper($_POST["GENin"]); // Retrieve the gene sequence
$TOTseq = $PROseq.$GENseq; // Find the coding strand we can target
$Template = complement($TOTseq); // Find the template strand we can target

echo "<div>
	The PAM sequence is $PAMseq.<br>
	The promoter sequence is $PROseq.<br>
	The gene sequence is $GENseq.<br>
	The coding strand is $TOTseq.<br>
	The template strand is $Template.<br><br>
</div>";

// NOT YET - create_db();

$len = 20; // The length of sgRNA to be formed
$res = willfind($Template, $PAMseq, $len); // Will we find results?

if ($res === true) { // If we find results, print them in a table
echo "<table>
	<th>Position of PAM</th>
	<th>Strand</th>
	<th>Following sequence to target</th>";

findsites($Template, $PAMseq, $len);

echo "</table>";
}
else { // If we don't find results, cry
	echo "No results found :(";
}

?>

</body>
</html>
