function y = complex2real(x)
    % Reinterprets complex/real array as real array
    %
    % SYNTAX
    %
    %   y = complex2real(x)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   x ... Full complex or real array
    %   y ... Full real array with interleaved real and imaginary entries of x
    %
    % NOTES
    %
    %   This is a less efficient implementation than the MEX implementation
    %   below, which is commented out as it might be broken
    assert(isnumeric(x) && ~issparse(x));
    if isreal(x)
        y = x;
    else
        shape = size(x);
        y = zeros([2*shape(1), shape(2:end)], class(x));
        y(1:2:end) = real(x);
        y(2:2:end) = imag(x);
    end
end


% FIXME: The efficient implementation (commented out):
%        SHAREDCHILD's complex2real() might be broken;
%        perhaps R2020a release of SHAREDCHILD will deliver a fix
function y = complex2real_(x)  %#ok<DEFNU>
    % Reinterprets complex/real array as real array (no new memory allocated)
    %
    % SYNTAX
    %
    %   y = complex2real(x)
    %
    % INPUT/OUTPUT PARAMETERS
    %
    %   x ... Full complex or real array
    %   y ... Full real array with interleaved real and imaginary entries of x
    %
    % NOTES
    %
    %   For real X, COMPLEX2REAL(X) returns X.
    %
    %   For complex X, COMPLEX2REAL(X) returns real Y with twice
    %   as many entries along the first dimension of X. Original
    %   real and imaginary entries of X(:) are interleaved in Y(:).
    %
    %   No copy is created, data are shared.
    %
    %   Roughly equivalent to the following code (except that no
    %   copy is made):
    %
    %       function y = complex2real(x)
    %           assert(isnumeric(x) && ~issparse(x));
    %           if isreal(x)
    %               y = x;
    %           else
    %               shape = size(x);
    %               y = zeros([2*shape(1), shape(2:end)], class(x));
    %               y(1:2:end) = real(x);
    %               y(2:2:end) = imag(x);
    %           end
    %       end


    try
        % Call MEX implementation
        y = util.complex2real_(x);
    catch ME
        switch ME.identifier
        case {'MATLAB:undefinedVarOrClass', 'MATLAB:UndefinedFunction'}
            % MEX binary does not exist; compile and try again
            compile_complex2real_();
            y = util.complex2real_(x);
        otherwise
            rethrow(ME)
        end
    end
end


function compile_complex2real_()
    path = fileparts(mfilename('fullpath'));
    cfile = fullfile(path, 'complex2real_src', 'complex2real_.c');
    w = warning('off', 'MATLAB:mex:GccVersion');  % Mute warning
    mex(cfile, '-outdir', path, '-silent', '-R2018a');
    warning(w);  % Restore warning
end
