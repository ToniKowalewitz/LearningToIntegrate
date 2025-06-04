function ptr = getPr(var)
    % Returns memory address to real part of the given array.
    %
    % SYNTAX
    %
    %   ptr = getPr(var)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   var ... Array of type mxDOUBLE_CLASS.
    %   ptr ... Scalar of class 'uint32' or 'uint64' (depending on the
    %           architecture's pointer data type size) with the starting
    %           address of the real data of the array.
    %
    % REMARKS
    %
    %   See https://mathworks.com/help/matlab/apiref/mxgetpr.html.

    try
        % Call MEX implementation
        ptr = util.getPr_(var);
    catch ME
        switch ME.identifier
        case {'MATLAB:undefinedVarOrClass', 'MATLAB:UndefinedFunction'}
            % MEX binary does not exist; compile and try again
            compile_getPr_();
            ptr = util.getPr_(var);
        otherwise
            rethrow(ME)
        end
    end
end


function compile_getPr_()
    path = fileparts(mfilename('fullpath'));
    cfile = fullfile(path, 'getPr_.c');
    w = warning('off', 'MATLAB:mex:GccVersion');  % Mute warning
    mex(cfile, '-outdir', path, '-silent');
    warning(w);  % Restore warning
end
