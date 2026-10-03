function varargout = run(name, varargin)
%RUN The only way to call a fabric estimator.
%
% varargout = ptt.estimators.run(name, ...)
%
% Every fabric estimator lives in +ptt/+estimators/private/, and MATLAB lets
% nothing outside this folder call a function there: not ptt.<name>, not a
% bare <name>, not str2func('<name>'). This function is the one door. It
% runs exactly ONE estimator, the one named, with the arguments given, and
% it refuses unless that estimator belongs to the method selected for this
% session with ptt.estimators.use. Two methods therefore cannot be mixed
% into one product by accident: switching method is a visible
% ptt.estimators.use call in the code that does it.
%
% The first scalar struct output is stamped with a field .estimator holding
% the method, the estimator's name and its convention - which axis an angle
% names, in which frame, and whether dlam is signed - taken from
% ptt.estimators.registry. Products built from the output carry that stamp
% forward, so what produced a number never has to be inferred.
%
% Inputs
%   name      an estimator name listed by ptt.estimators.registry
%   varargin  that estimator's own arguments, unchanged
%
% Example
%   ptt.estimators.use('quadpol_ls');
%   out = ptt.estimators.run('quadpolFabricLS', S, z, opts);
%
% See also ptt.estimators.use, ptt.estimators.current,
%   ptt.estimators.registry.

if nargin < 1 || ~((ischar(name) && isrow(name)) || (isstring(name) && isscalar(name)))
  error('ptt:estimators:name', ...
    'name the estimator: ptt.estimators.run(name, ...)');
end
name = char(name);
[method, convention] = H_lookup(name);
selected = selectedMethod();
if isempty(selected)
  error('ptt:estimators:noMethod', ...
    ['no estimation method is selected. Call ptt.estimators.use(''%s'') ' ...
     'first: %s belongs to that method.'], method, name);
end
if ~strcmp(selected, method)
  error('ptt:estimators:notInMethod', ...
    ['%s belongs to method ''%s'', but this session selected ''%s''. ' ...
     'One method per run: call ptt.estimators.use(''%s'') to switch ' ...
     'deliberately.'], name, method, selected, method);
end

impl = str2func(name);      % resolved here, where private/ is visible
[varargout{1:max(nargout, 1)}] = impl(varargin{:});

stamp = struct('method', method, 'name', name, 'convention', convention);
for k = 1:numel(varargout)
  if isstruct(varargout{k}) && isscalar(varargout{k})
    varargout{k}.estimator = stamp;
    break
  end
end
end

function [method, convention] = H_lookup(name)
R = ptt.estimators.registry();
for i = 1:numel(R)
  hit = strcmp({R(i).members.name}, name);
  if any(hit)
    method = R(i).method;
    convention = R(i).members(hit).convention;
    return
  end
end
names = arrayfun(@(r) strjoin({r.members.name}, ', '), R, 'UniformOutput', false);
error('ptt:estimators:unknown', ...
  'no estimator named ''%s''. The estimators are: %s', name, strjoin(names, '; '));
end
