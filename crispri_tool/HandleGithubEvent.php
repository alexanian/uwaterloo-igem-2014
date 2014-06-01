$JSONString = file_get_contents('php://input');
$JSONObject = json_decode($input);

if( json_last_error() === JSON_ERROR_NONE )
{
    if( fnmatch( "*/tool_beta", $JSONString['ref'], FNM_PATHNAME ) )
    {
        $GitResult = `git pull -f`;
    }
}
else
{
    header('HTTP/1.0 404 Not Found');
}
