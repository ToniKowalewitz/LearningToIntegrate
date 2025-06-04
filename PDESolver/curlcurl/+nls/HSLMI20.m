classdef HSLMI20 < solving.HSLMI20

    methods (Access = public)

        function obj = HSLMI20(varargin)
            warning('HSLMI20:deprecation', ...
                    'nls.HSLMI20 is deprecated. Use solving.HSLMI20 instead.');
            obj = obj@solving.HSLMI20(varargin{:});
        end

        function [x, inform] = precondition(obj, b)
            % Solve A*x = b
            %
            % This is a deprecated method supporting only single
            % vector b, which can be a column or a row. Replaced
            % by HSLMI20.solve(), which is vectorized.

            warning('HSLMI20:deprecation', ...
                    'Method HSLMI20.precondition() is deprecated. Use solve() instead.');
            if ~( iscolumn(b) || isrow(b) )
                error('Expected a vector');
            end
            [x, inform] = obj.solve(b(:));
        end

    end

end
