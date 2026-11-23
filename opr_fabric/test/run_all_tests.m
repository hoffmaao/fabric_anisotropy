function varargout = run_all_tests(selection)
%RUN_ALL_TESTS Run the test suite; error if any test fails.
%
%   run_all_tests              every test_*.m in this directory
%   run_all_tests('quick')     all but the slow ones - the smoke test for a
%                              fresh clone, and what the gate runs
%   run_all_tests({'test_node_fabric','test_xcorr_traveltime'})   just those
%   results = run_all_tests(...)   also return name/status/seconds/message
%
% From a shell, where any failure becomes a non-zero exit status:
%
%   matlab -batch "addpath('<clone>/opr_fabric/test'); run_all_tests quick"
%
% or run_all_tests.sh beside this file, which finds MATLAB for you.
%
% Each test is a script that starts with clear and ends by erroring if any
% of its checks failed, so a thrown error is the whole verdict. A test runs
% inside a helper function's workspace, where its clear cannot reach the
% runner or the caller's base workspace. After it, the path, working
% directory and warning state are restored and the figures it opened are
% closed (figures that were already open stay), so nothing leaks into the
% next test. A failing test is reported and the suite carries on, so one
% run names every failure.

here = fileparts(mfilename('fullpath'));
listing = dir(fullfile(here, 'test_*.m'));
names = sort(erase({listing.name}, '.m'));

% Tests that take over two minutes each on a laptop go here, and 'quick'
% leaves them out. There are none now: the three tests take a little over
% 2 min together on an Apple-silicon laptop (R2025b; test_node_fabric
% 1.7 min, test_xcorr_traveltime 0.5 min) and about 10 min on mem1 at 8
% threads. Update when a test is added or its cost changes, so 'quick'
% stays quick.
slow = {};

if nargin == 0 || isequal(selection, 'all')
  % every test
elseif isequal(selection, 'quick')
  names = setdiff(names, slow, 'stable');
else
  want = erase(cellstr(selection), '.m');
  unknown = setdiff(want, names);
  if ~isempty(unknown)
    error('run_all_tests:unknown', 'No such test in %s: %s', here, ...
      strjoin(unknown, ', '));
  end
  names = want;
end

results = struct('name', names, 'status', '', 'seconds', NaN, 'message', '');
t_all = tic;
for k = 1:numel(names)
  fprintf('\n==== %s (%d of %d) ====\n', names{k}, k, numel(names));
  saved_path = path;
  saved_dir = pwd;
  saved_warn = warning;
  saved_figs = findall(groot, 'Type', 'figure');
  t = tic;
  try
    H_run(fullfile(here, [names{k} '.m']));
    results(k).status = 'PASS';
  catch err
    results(k).status = 'FAIL';
    results(k).message = err.message;
    fprintf(2, '%s\n', getReport(err, 'extended', 'hyperlinks', 'off'));
  end
  results(k).seconds = toc(t);
  path(saved_path);
  cd(saved_dir);
  warning(saved_warn);
  figs = findall(groot, 'Type', 'figure');
  close(figs(~ismember(figs, saved_figs)), 'force');
end

fprintf('\n%-30s %-6s %9s\n', 'test', 'result', 'seconds');
for k = 1:numel(results)
  fprintf('%-30s %-6s %9.1f  %s\n', results(k).name, results(k).status, ...
    results(k).seconds, results(k).message);
end
failed = strcmp({results.status}, 'FAIL');
fprintf('%d passed, %d failed (%.1f min)\n', nnz(~failed), nnz(failed), ...
  toc(t_all)/60);

if nargout > 0
  varargout{1} = results;
end
if any(failed)
  error('run_all_tests:failed', '%d of %d test(s) failed: %s', ...
    nnz(failed), numel(results), strjoin(names(failed), ', '));
end
end

function H_run(script_fn)
% A function workspace of its own for the test script: run evaluates the
% script here, so its clear empties only this workspace.
run(script_fn);
end
