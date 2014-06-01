{
    "Result" : "<?php

$JSONString = file_get_contents('php://input');
$JSONObject = json_decode($JSONString);

if( json_last_error() === JSON_ERROR_NONE )
{
    if( isset( $JSONObject->ref ) )
    {
        if( fnmatch( "*/tool_beta", $JSONObject->ref, FNM_PATHNAME ) )
        {
            echo `git pull -f`;
        }
        else
        {
            echo "Wasn't expecting the input given.";
        }
    }
    else
    {
        echo "Wasn't expecting the input given.";
    }
}
else
{
    echo "JSON object was invalid.";
}

?>"
}
