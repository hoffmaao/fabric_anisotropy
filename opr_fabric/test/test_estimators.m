%TEST_ESTIMATORS The estimator folder admits one method at a time, through one door.
%
% Every fabric estimator lives in +ptt/+estimators/private/ and is reached
% only through ptt.estimators.run, which runs the one estimator named,
% refuses any estimator outside the method selected with
% ptt.estimators.use, and stamps what ran on its output. This exercises that
% contract from outside the folder, where every other caller sits.
%
%   1. NO WAY AROUND THE DOOR. An estimator cannot be called as ptt.<name>,
%      as ptt.estimators.<name>, bare, or through str2func.
%   2. ONE METHOD. With no method selected the door refuses; with one
%      selected it refuses an estimator that belongs to another; an
%      unknown method or estimator is refused by name.
%   3. THE STAMP. Through the door an estimator returns its usual outputs,
%      and its first struct output gains .estimator: the method, the name
%      and the convention the registry records - on the struct even when a
%      numeric output comes first.
%   4. THE REGISTRY IS THE FOLDER. Every estimator belongs to exactly one
%      method, every registered estimator resolves through the door (with
%      no arguments it fails on its own missing inputs, inside the private
%      folder, not in the door), and no file in the private folder is left
%      unregistered.
%
% Run: matlab -batch "run('opr_fabric/test/test_estimators.m')"
clear;
thisDir = fileparts(mfilename('fullpath'));
root = fullfile(thisDir, '..', '..');
addpath(root);
fails = 0;

%% 1. no way around the door
probes = {
  'ptt.quadpolFabricLS',            @() ptt.quadpolFabricLS(1);
  'ptt.estimators.quadpolFabricLS', @() ptt.estimators.quadpolFabricLS(1);
  'bare quadpolFabricLS',           @() quadpolFabricLS(1);
  'str2func(''quadpolFabricLS'')',  @() feval(str2func('quadpolFabricLS'), 1);
  };
for i = 1:size(probes, 1)
  id = H_err(probes{i, 2});
  ok = any(strcmp(id, {'MATLAB:undefinedVarOrClass', 'MATLAB:UndefinedFunction'}));
  fprintf('1. %-34s refused (%s): %s\n', probes{i, 1}, id, H_tick(ok));
  fails = fails + ~ok;
end

%% 2. one method at a time
ptt.estimators.use('none');
id = H_err(@() ptt.estimators.run('traveltimeFabricML', 1, 1));
ok = strcmp(id, 'ptt:estimators:noMethod');
fprintf('2. nothing selected, the door refuses (%s): %s\n', id, H_tick(ok));
fails = fails + ~ok;

prev = ptt.estimators.use('quadpol_ls');
ok = isempty(prev) && strcmp(ptt.estimators.current(), 'quadpol_ls');
fprintf('   use returns the previous selection and current reports the new: %s\n', H_tick(ok));
fails = fails + ~ok;
id = H_err(@() ptt.estimators.run('ershadiFabric', 1, 1));
ok = strcmp(id, 'ptt:estimators:notInMethod');
fprintf('   quadpol_ls selected, an ershadi estimator is refused (%s): %s\n', id, H_tick(ok));
fails = fails + ~ok;
id = H_err(@() ptt.estimators.use('best_one'));
ok = strcmp(id, 'ptt:estimators:unknownMethod');
fprintf('   an unknown method is refused (%s): %s\n', id, H_tick(ok));
fails = fails + ~ok;
id = H_err(@() ptt.estimators.run('selectedMethod'));
ok = strcmp(id, 'ptt:estimators:unknown');
fprintf('   the folder''s own state is not an estimator (%s): %s\n', id, H_tick(ok));
fails = fails + ~ok;

%% 3. the stamp
R = ptt.estimators.registry();
conv = @(m, n) R(strcmp({R.method}, m)).members(strcmp({R(strcmp({R.method}, m)).members.name}, n)).convention;
ptt.estimators.use('nymand');
z = (30:10:600).';
dtau = 0.0637 * cumsum(0.02 + 0.03 * (z - 30) / 570) * 10;
o = ptt.estimators.run('traveltimeFabricML', dtau, z, struct('prior', [0.03 0.05 100], 'sigma_ns', 0.03));
ok = isfield(o, 'estimator') && strcmp(o.estimator.method, 'nymand') ...
  && strcmp(o.estimator.name, 'traveltimeFabricML') ...
  && strcmp(o.estimator.convention, conv('nymand', 'traveltimeFabricML')) ...
  && isfield(o, 'dlam') && numel(o.dlam) == numel(z);
fprintf('3. one output: the usual fields plus .estimator (nymand / traveltimeFabricML): %s\n', H_tick(ok));
fails = fails + ~ok;

ptt.estimators.use('copol');
par = ptt.defaultParams();
par.lam_x_sfc = 1/3; par.lam_z_sfc = 1/3; par.lam_x_bed = 0.05; par.lam_z_bed = 0.50;
obs = struct('L', 5, 'z', ([0.9 0.8 0.65 0.5 0.35 0.2] * par.H).');
obs.dtau = ptt.twttDifference(par, obs.L * ones(size(obs.z)), obs.z);
[dl, out] = ptt.estimators.run('invertHorizontalFabric', obs, par);
ok = isnumeric(dl) && numel(dl) == numel(obs.z) && isfield(out, 'estimator') ...
  && strcmp(out.estimator.name, 'invertHorizontalFabric');
fprintf('   two outputs: dlam stays numeric, the struct after it is stamped: %s\n', H_tick(ok));
fails = fails + ~ok;

%% 4. the registry is the folder
names = arrayfun(@(r) {r.members.name}, R, 'UniformOutput', false);
names = [names{:}];
ok = numel(unique(names)) == numel(names);
fprintf('4. %d estimators in %d methods, none in two: %s\n', numel(names), numel(R), H_tick(ok));
fails = fails + ~ok;
bad = {};
privdir = fullfile('+ptt', '+estimators', 'private');
for i = 1:numel(R)
  ptt.estimators.use(R(i).method);
  for n = {R(i).members.name}
    % Called with no arguments, a resolved estimator fails on its OWN
    % missing inputs, so the error is raised inside the private folder; an
    % unresolved name would fail in the door, at the call to it.
    [id, where] = H_err(@() ptt.estimators.run(n{1}));
    if isempty(id) || ~contains(where, privdir)
      bad{end+1} = sprintf('%s (%s at %s)', n{1}, id, where); %#ok<SAGROW>
    end
  end
end
ok = isempty(bad);
fprintf('   every registered estimator resolves through the door: %s%s\n', H_tick(ok), ...
  H_list(bad));
fails = fails + ~ok;
files = dir(fullfile(root, '+ptt', '+estimators', 'private', '*.m'));
files = setdiff(erase({files.name}, '.m'), {'selectedMethod'});
orphans = setdiff(files, names);
missing = setdiff(names, files);
ok = isempty(orphans) && isempty(missing);
fprintf('   the private folder holds exactly the registered estimators: %s%s%s\n', ...
  H_tick(ok), H_list(orphans), H_list(missing));
fails = fails + ~ok;

ptt.estimators.use('none');
fprintf('\n%s\n', H_tick(fails == 0));
if fails > 0
  error('test_estimators:failed', '%d check(s) failed', fails);
end

function [id, where] = H_err(f)
% the error identifier f raises and the file it was raised in ('' if none)
id = ''; where = '';
try
  f();
catch e
  id = e.identifier;
  if ~isempty(e.stack), where = e.stack(1).file; end
end
end

function s = H_list(c)
if isempty(c), s = ''; else, s = sprintf(' [%s]', strjoin(c, ', ')); end
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
