$.ajaxSetup ({
    cache: false
});

$.fn.serializeObject = function()
{
    var o = {};
    var a = this.serializeArray();
    $.each(a, function() {
        if (o[this.name] !== undefined) {
            if (!o[this.name].push) {
                o[this.name] = [o[this.name]];
            }
            o[this.name].push(this.value || '');
        } else {
            o[this.name] = this.value || '';
        }
    });
    return o;
};

$(document).ready(function() {
	// Hide results container divs on start-up
	$('#resultsContainer').hide();
	$('#plasmidContainer').hide();
	$('#results').dataTable();
	drawArc();
	drawCircle();	
	
	// What to do on submitting the search button
	$("#inputData").submit(function(event) {
		event.preventDefault();
			
		// If a custom PAM sequence has been entered, set the value of
		// the custom radio button to the input
		if ( $('#customPAM').val() ) {
			$('#pamseq3').val($('#customPAM').val().toUpperCase());
		}
		
		$('#results').dataTable().fnDestroy(); 
		$("#results").dataTable( {
		    "processing" : true,
		    "serverSide" : true,
		    "ajax" : {
			    type: "POST",
			    url: this.action,
			    data: $("#inputData").serializeObject()
		        },
	        "columns": [
                    { "data" : "GeneName" },
                    { "data" : "Strand" },
                    { "data" : "Position" },
                    { "data" : "sgRNA" }
                ]
	        } );
		
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
