function previous = use(method)
%USE Select the ONE estimation method this session may call.
%
% previous = ptt.estimators.use(method)
%
% ptt.estimators.run refuses any estimator outside the selected method. A
% script, test or pipeline selects its method once, before its first call;
% one that genuinely compares two methods switches with a second call, so
% the switch is visible where it happens. Returns the method that was
% selected before ('' if none), so a caller can put it back.
%
% use('none') deselects, after which ptt.estimators.run refuses
% everything; the test runner does this before every test, so a test that
% forgets to select its own method fails instead of inheriting the last
% test's.
%
% The selection lives for the MATLAB session. `clear functions` and
% `clear all` reset it, after which ptt.estimators.run errors until a
% method is selected again - a reset can make a call fail, never silently
% change which method answers.
%
% See also ptt.estimators.run, ptt.estimators.current,
%   ptt.estimators.registry.

if nargin < 1 || ~((ischar(method) && isrow(method)) || (isstring(method) && isscalar(method)))
  error('ptt:estimators:method', ...
    'name the method: ptt.estimators.use(method)');
end
method = char(method);
if strcmp(method, 'none')
  previous = selectedMethod();
  selectedMethod('');
  return
end
R = ptt.estimators.registry();
if ~any(strcmp({R.method}, method))
  error('ptt:estimators:unknownMethod', ...
    'no method named ''%s''. The methods are: %s', method, ...
    strjoin({R.method}, ', '));
end
previous = selectedMethod();
selectedMethod(method);
end
