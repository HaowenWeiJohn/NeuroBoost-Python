function env = loadenv(envFile)
%LOADENV  Parse a .env file and set key/value pairs with setenv.
%
%   ENV = LOADENV(FILE)               – explicit .env path (preferred)
%   ENV = LOADENV()                   – looks for ".env" in pwd
%
%   If FILE does *not* exist, LOADENV returns an empty struct and
%   leaves the workspace untouched (no error).

    if nargin == 0 || isempty(envFile)
        envFile = '.env';
    end

    if ~isfile(envFile)                 % ← changed behaviour
        warning('loadenv:NotFound', ...
                '.env file "%s" not found – skipped.', envFile);
        env = struct();                 % nothing loaded
        return
    end

    txt  = fileread(envFile);
    rows = regexp(txt, '\r?\n', 'split');

    env = struct();
    for r = string(rows)
        r = strip(r);
        if r=="" || startsWith(r,"#"), continue, end
        tokens = regexp(r, '([^=]+)=(.*)', 'tokens', 'once');
        if isempty(tokens),  continue,  end
        key = strtrim(tokens{1});
        val = regexprep(strtrim(tokens{2}), '^["'']?(.*?)["'']?$', '$1');
        setenv(key, val);
        env.(matlab.lang.makeValidName(key)) = val;
    end
end
