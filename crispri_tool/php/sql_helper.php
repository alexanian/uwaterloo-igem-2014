<?php

//---------------------------------------------------
function connect_db() {
	
$user = 'root';
$pass = '';
$db = 'possible_targets';

$mysqli = new mysqli('localhost',$user,$pass,$db);

if (mysqli_connect_error()) {
    die('Connect Error (' . mysqli_connect_errno() . ') ' . mysqli_connect_error());
}

// echo 'Successful connection...' . $mysqli->host_info . "<br><br>";

$trunc = "TRUNCATE TABLE Targets";
if (!mysqli_query($mysqli,$trunc)) {
	die('Error: ' . mysqli_error($mysqli));
}

// Create table
$sql="CREATE TABLE Targets(Gene VARCHAR(10), Strand VARCHAR(10), Position INT, Result TINYTEXT, Repression DOUBLE(20,5))";
mysqli_query($mysqli,$sql);

}

//---------------------------------------------------
function insert_row($gene, $strand, $pos, $sgRNA) {

$user = 'root';
$pass = '';
$db = 'possible_targets';

$conn = mysqli_connect('localhost',$user,$pass,$db);
$new_record = "INSERT INTO Targets (Gene, Strand, Position, Result, Repression)
	VALUES ('$gene', '$strand', '$pos', '$sgRNA', '0')";

if (!mysqli_query($conn,$new_record)) {
  die('Error: ' . mysqli_error($conn));
}

}


?>