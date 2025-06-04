%#ok<*INUSD>

function tests = test_style
  % Run all function tests in this file
  tests = functiontests(localfunctions);
end


function test_run_linter(t)
    files = dir('+*/**/*.m');

    for f=files'
        fn = strcat(f.folder, '/', f.name);

        results = mlint(fn, '-id');  %#ok<*MLNT>

        if ~isempty(results)
            t.verifyFail(strcat('Failed linting:', 32, fn));
            for r=results'
                fprintf('L %d (C %d-%d): %s: %s\n', r.line, r.column(:), r.id, r.message);
            end
        end
    end
end


function test_whitespace(t)
    cmd = 'git diff-index --check $(git hash-object -t tree /dev/null)';
    util.run_sys_cmd(cmd);
end


function test_codestyle(t)
end
