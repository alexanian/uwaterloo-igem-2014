{
    "Result" : "<?php

$JSONString = file_get_contents('php://input');
$JSONObject = json_decode($JSONString);

if( json_last_error() === JSON_ERROR_NONE )
{
    if( isset( $JSONObject->ref ) )
    {
        if( fnmatch( "*/tool_beta", $JSONObject->ref ) )
        {
            echo `git pull -f`;
        }
        else
        {
            echo "Not going to publish, wasn't a push to tool_beta branch.";
        }
    }
    else
    {
        echo "Wasn't expecting the input given. Missing ref.";
    }
}
else
{
    echo "JSON object was invalid.";
}

?>"
}
