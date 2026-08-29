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
% username into 22 MATLAB scripts and 6 shell launchers, which is what made
% this code unrunnable by anyone else on the CReSIS machines.
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
% exactly the override the variable exists for.
%
% FABRIC_ROOT is therefore checked here for exactly what the derived value
% gives for free: it must name a directory that EXISTS and must resolve to
% an ABSOLUTE path, a derived work_root being both by construction (it is
% the parent of the directory this file was found in), and it is rejected
% rather than created, since a typo created on first write would quietly
% collect the products. Why an override is held to the derived value's test
% rather than a weaker one is an argument about both helpers, so it is
% owned by docs/scripts.md, "Running on the CReSIS machines".
%
% The part of that check which is MATLAB's alone: the override is resolved
% by ENTERING it, not via what() - what() consults the MATLAB search path,
% so a relative override could bind to a directory other than the one
% isfolder just validated relative to the working directory, which is the
% whole failure this check exists to stop.
%
% fabric_paths.sh beside this file is the shell launchers' mirror, anchored
% on the same +ptt walk. It diverges deliberately in two ways, and what
% each means on THIS side is: stages/ is not required here, only there; and
% the work_root handed back here has its symlinks RESOLVED, because MATLAB's
% pwd reports the OS getcwd, where bash's logical pwd returns the path as
% written. Both divergences are claims about the pair, so docs/scripts.md,
% "Running on the CReSIS machines", owns the reasoning for each and neither
% header restates it.
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
  given = work_root;
  if ~isfolder(given)
    error('fabric_paths:noRoot', ...
      ['FABRIC_ROOT is %s, which is not a directory - a root that does ' ...
      'not exist would put every product in a tree nothing else reads'], ...
      given);
  end
  try
    here = cd(given);
  catch err
    error('fabric_paths:noRoot', ...
      ['FABRIC_ROOT is %s, which cannot be entered - check permissions ' ...
      '(%s)'], given, err.message);
  end
  restore = onCleanup(@() cd(here));
  work_root = pwd;
  clear restore;
  if isempty(work_root)
    error('fabric_paths:noRoot', ...
      'FABRIC_ROOT is %s, which did not resolve to an absolute path', given);
  end
end

end
