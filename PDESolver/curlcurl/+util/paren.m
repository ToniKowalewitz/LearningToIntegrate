function out = paren(x, varargin)
    % Returns a cell including the element 'varargin' of array 'x'.
    %
    % SYNTAX
    %
    %   out = paren(x, varargin)
    %
    % INPUT PARAMETERS
    %   x ... Array.
    %
    % OUTPUT PARAMETERS
    %   out      ... Selected element.
    %   varargin ... Index of the argument that should be returned. E.g.
    %                specify '2' to assign {x{2}}.

    out = x(varargin{:});
end
