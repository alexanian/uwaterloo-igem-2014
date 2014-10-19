// Draw plasmid
function drawArc(){
	var canvas = document.getElementById('cartoonPlasmid');
	var context = canvas.getContext('2d');
	var x = canvas.width / 2;
	var y = canvas.height / 2;
	var radius = 75;
	var startAngle = -0.35 * Math.PI;
	var endAngle = 0.35 * Math.PI;

	context.beginPath();
	context.arc(x, y, radius, startAngle, endAngle);
	context.lineWidth = 15;
	
	context.strokeStyle = 'blue';
	context.stroke();
}

// Draw circle
function drawCircle(){
	
	var canvas = document.getElementById('cartoonPlasmid');
	var context = canvas.getContext('2d');
	var centerX = canvas.width / 2;
	var centerY = canvas.height / 2;
	var radius = 75;
	var startAngle = 0.35 * Math.PI;
	var endAngle = -0.35 * Math.PI;
	
	context.beginPath();
	context.arc(centerX, centerY, radius, startAngle, endAngle);
	//context.fillStyle = 'green';
	//context.fill();
	context.lineWidth = 8;
	
	context.strokeStyle = 'black';
	context.stroke();
	
	
}


