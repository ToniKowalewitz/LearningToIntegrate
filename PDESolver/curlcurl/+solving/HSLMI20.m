classdef HSLMI20 < handle

    properties (GetAccess = public, SetAccess = private)
        inform
        handle
    end

    methods (Access = public)

        function obj = HSLMI20(A, control)
            % Construct HSLMI20 AMG solver
            %
            % INPUT ARGS
            %   A ... sparse square matrix
            % OPTIONAL ARGS
            %   control ... HSLMI20 parameters obtained with hsl_mi20_control()

            if ~obj.is_installed()
                error('HSL_MI20 package not installed!');
            end

            hsl_mi20_startup();

            if nargin == 1
                control = hsl_mi20_control();
            end

            [obj.inform, obj.handle] = hsl_mi20('setup', A, control);
        end

        function delete(obj)
            fprintf('Deleting hsl_mi20 handle %d\n', obj.handle);
            hsl_mi20('destroy', obj.handle);
        end

        function [x, inform] = solve(obj, b)
            % Solve equation A*x = b
            %
            % b can have multiple columns

            % Allocate return values
            x = zeros(size(b));
            inform = obj.create_inform_struct(size(b, 2));

            % Solve column-by-column
            for j = 1:size(b, 2)
                [x(:, j), inform(j)] = hsl_mi20('precondition', b(:, j), obj.handle);
            end
        end

    end

    methods (Static)

        function result = is_installed()
            try
                hsl_mi20_startup();
                result = true;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                    result = false;
                else
                    rethrow(ME);
                end
            end
        end

        function inform = create_inform_struct(n)
            inform(1:n) = struct('flag', -Inf, 'clevels', 0, 'cpoints', 0, 'cnnz', 0);
        end

    end

end
