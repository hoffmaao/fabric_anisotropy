function [code_root, work_root] = fabric_paths()
%FABRIC_PATHS Repo root and scratch work root, derived from this file.
%
% [code_root, work_root] = fabric_paths()
%
% Processing happens in a scratch work root laid out as
%
%   <work>/code      this repo (code_root)
%   <work>/stages    its products
%   <work>/invert_logs, locks, ...
%
% and every server script needs those two directories. They are DERIVED
% from this file's own location rather than named, so the same tree runs
% from any user's scratch and from a second checkout used to test a change
% without disturbing the live work root. Hardcoding them baked one
% username into 21 scripts, which is what made this code unrunnable by
% anyone else on the CReSIS machines.
%
% Deliberately not relative to the WORKING DIRECTORY. The batch launchers
% cd to <work>, but the one-liner form
%
%   matlab -batch "addpath('<code>/opr_fabric/server'); \
%                  day_seg='20250108_02'; frm=9; run_season_ridge_a"
%
% needs this directory on the path already, to find the script at all - so
% the cwd that lets MATLAB find a script is not necessarily the cwd that
% would make './code' resolve. This file's location is fixed by the deploy
% layout and is the same under both.
%
% There is no bootstrap problem in putting this in a function: it sits
% beside the scripts that call it, so whatever put them within MATLAB's
% reach put this here too.
%
% DATA roots stay absolute and are NOT derived here - /cresis/... are real
% mount points on the CReSIS machines, not part of the work tree, and each
% season declares its own site_root. So does shared software such as
% /kucresis/scratch/software/snaphu.
%
% The walk looks for the +ptt toolbox rather than counting directories, so
% it stays right if the tree is nested differently, and it fails loudly
% rather than putting a wrong directory on the path and dying later on a
% missing function. FABRIC_ROOT overrides work_root - and ONLY work_root -
% for a layout whose products do not sit beside the code. code_root always
% comes from this file, so callers must take it from here rather than
% rebuilding it as fullfile(work_root, 'code'), which is wrong under
% exactly the override the variable exists for. The override must name a
% directory that exists, and is resolved to an absolute path: a typo would
% otherwise be created on first write and quietly collect the products.
%
% fabric_paths.sh beside this file is the shell launchers' mirror of the
% same two rules, anchored on the same +ptt walk. It differs in one stated
% way, which its header gives: it requires stages/ to already exist.
%
% See also run_quadpol_pipeline, run_season_frame.

d = fileparts(mfilename('fullpath'));
code_root = '';
for k = 1:5
  if exist(fullfile(d, '+ptt'), 'dir') == 7
    code_root = d;
    break
  end
  parent = fileparts(d);
  if isempty(parent) || strcmp(parent, d)
    break
  end
  d = parent;
end
if isempty(code_root)
  error('fabric_paths:noToolbox', ...
    ['no +ptt toolbox above %s - the server scripts must sit inside the ' ...
    'repo tree (<work>/code) so the toolbox they need can be found'], ...
    mfilename('fullpath'));
end

work_root = getenv('FABRIC_ROOT');
if isempty(work_root)
  work_root = fileparts(code_root);
else
  if ~isfolder(work_root)
    error('fabric_paths:noRoot', ...
      ['FABRIC_ROOT names %s, which is not a directory - an override ' ...
      'that does not exist would put every product in a tree nothing ' ...
      'else reads'], work_root);
  end
  w = what(work_root);
  work_root = w(1).path;
end

end
