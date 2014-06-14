$.ajaxSetup ({
    cache: false
});

$(document).ready(function() {
	// What to hide on start-up
	$('#resultsContainer').hide();
	$('#plasmidContainer').hide();
	drawArc();
	drawCircle();
	
	
	// What to do on submit
	$("#inputData").submit(function(event) {
		event.preventDefault();
		
		$.ajax({
			type: "POST",
			url: this.action,
			data: $("#inputData").serialize(), //serializes the form's elements.
			success: function(data){
				//alert(data); // show response from the php script.
				$("#phpHere").html(data);
			}
		});
		
		$("#inputContainer").hide();
		$("#resultsContainer").show();
		
		return false; // prevent to execute the actual submit of the form.
	});
	
	
	// Go back to entry
	$("#returnInput").click(function(){
		$("#inputContainer").show();
		$("#resultsContainer").hide();
	});
	
	
	// See plasmid view
	$("#viewPlasmid").click(function(){
		$("#plasmidContainer").show();
		$("#resultsContainer").hide();
	});
	
	
	// See plasmid view
	$("#returnTable").click(function(){
		$("#resultsContainer").show();
		$("#plasmidContainer").hide();
	});
	
});