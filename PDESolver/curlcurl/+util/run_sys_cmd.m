function out = run_sys_cmd(cmd, ignore)
    % Run system command with error handling.

    if nargin < 2
        ignore = false;
    end

    [ierr, out] = system(cmd);
    msg_text = 'command "%s" failed with return code %d and message:\n%s';
    if ierr > 0
        if ~ignore
            error(msg_text, cmd, ierr, out);
        else
            warning(msg_text, cmd, ierr, out);
        end
    end
end
