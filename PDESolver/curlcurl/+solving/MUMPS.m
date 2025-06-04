classdef MUMPS < handle

    properties (Access = public)
        id
    end

    properties (GetAccess = public, SetAccess = private)
        mumps
        A
    end

    methods (Access = public)

        function obj = MUMPS(A, id)
            % Factor sparse square matrix using MUMPS
            %
            % INPUT ARGS
            %   A ... sparse square matrix
            % OPTIONAL ARGS
            %   id ... MUMPS parameters obtained with initmumps()

            if ~obj.is_installed()
                error('MUMPS package not installed!');
            end

            obj.A = A;

            if isreal(obj.A)
                obj.mumps = @dmumps;
            else
                obj.mumps = @zmumps;
            end

            % Update MUMPS control structure
            obj.id = initmumps();
            if nargin > 1
                obj.id = update_struct(obj.id, id);
            end

            % Initialize MUMPS instance
            obj.id.JOB = -1;
            obj.id = obj.mumps(obj.id);
            obj.check_err();

            % Factorize (symbolic and numeric)
            obj.id.JOB = 4;
            obj.id = obj.mumps(obj.id, obj.A);
            obj.check_err();

            % Clean up
            obj.id.PIVNUL_LIST = -9999;
            obj.id.SYM_PERM = -9999;
            obj.id.UNS_PERM = -9999;
            obj.id.ROWSCA = -9999;
            obj.id.COLSCA = -9999;
            obj.id.SCHUR = -9999;
            obj.id.VAR_SCHUR = -9999;
        end

        function delete(obj)
            fprintf('Deleting MUMPS instance %d\n', obj.id.INST);
            obj.id.JOB = -2;
            obj.id = obj.mumps(obj.id);
            obj.check_err();
        end

        function [x, id] = solve(obj, b)
            % Solve equation A*x = b
            %
            % b can have multiple columns

            % Solve
            obj.id.RHS = b;
            obj.id.JOB = 3;
            obj.id = obj.mumps(obj.id, obj.A);
            obj.check_err();

            % Assign return values
            x = obj.id.SOL;
            if nargout > 1
                id = obj.id;
            end

            % Clean up
            obj.id.SOL = -9999;
            obj.id.RHS = -9999;
            obj.id.REDRHS = -9999;
        end

        function check_err(obj)
            if isempty(obj.id) || obj.id.INFOG(1) == 0
                return
            end

            msg = sprintf('INFOG(1) = %d\nINFOG = %s', ...
                          obj.id.INFOG(1), mat2str(obj.id.INFOG));

            if obj.id.INFOG(1) > 0
                warning(['MUMPS WARNING: ', msg]);
            elseif obj.id.INFOG(1) < 0
                error(['MUMPS ERROR: ', msg]);
            end
        end

    end

    methods (Static)

        function result = is_installed()
            try
                [~] = initmumps();
                result = true;
            catch ME
                if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
                    result = false;
                else
                    rethrow(ME);
                end
            end
        end

    end

end


function s = update_struct(s, s_new)

    % Remove fields in s_new from s,
    % raise error if field not in s
    s = rmfield(s, fieldnames(s_new));

    % Concatenate s and s_new
    fields = [fieldnames(s); fieldnames(s_new)];
    s = cell2struct([struct2cell(s); struct2cell(s_new)], fields, 1);
end
