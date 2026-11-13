% run_tests: run the four self-tests (lake at rest, conservation, AMR vs uniform, 2D); returns true if all pass
function ok = run_tests()
    setup_paths();
    names = {'test_lib_a','test_lib_b','test_amr','test_swe2d'};
    results = false(size(names));
    times = zeros(size(names));

    for i = 1:numel(names)
        try
            tstart = tic;
            run_one(names{i});
            times(i) = toc(tstart);
            results(i) = true;
        catch ME
            times(i) = toc(tstart);
            fprintf('Error in %s:\n%s\n', names{i}, ME.message);
            if ~isempty(ME.stack)
                frame = ME.stack(1);
                fprintf('  at %s (line %d)\n', frame.name, frame.line);
            end
        end
    end

    fprintf('Test results:\n');
    fprintf('%-12s %-5s %s\n', 'Name', 'Pass', 'Seconds');
    for i = 1:numel(names)
        fprintf('%-12s %-5s %.4f\n', names{i}, ternary(results(i),'PASS','FAIL'), times(i));
    end

    ok = all(results);
end

function run_one(name)
    run(name);
end

function out = ternary(c,t,f)
    if c, out = t; else, out = f; end
end
