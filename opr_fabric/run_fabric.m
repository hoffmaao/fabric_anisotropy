% script run_fabric
%
% Script for running fabric.m (usually just used for debugging).
%
% Runs straight from a fresh clone of this repo, wherever it sits (e.g.
% /kucresis/scratch/<user>/scripts/fabric_anisotropy) and whoever owns it:
% the Paths section below puts this checkout's +ptt and fabric.m /
% fabric_task.m on the MATLAB path, so nothing is copied into the OPR tree
% or run_opr. The only prerequisite is OPR itself, set up as usual by your
% startup.m (gRadar). Then, from any directory:
%
%   run('<clone>/opr_fabric/run_fabric.m')
%
% Requires CSARP_polarimetric products from polarimetric.m (run with
% coregistration.en = true and ideally snaphu_en = true); fabric.in_path
% below names which one.
%
% See also: run_master.m, master.m, fabric.m, fabric_task.m,
% run_polarimetric.m, polarimetric.m

%% Paths (automatic)
% =====================================================================
% This checkout's root holds +ptt and this directory holds fabric.m and
% fabric_task.m; both go to the FRONT of the path so an older copy
% elsewhere (e.g. one deployed into opr/matlab or run_opr) cannot shadow
% them. Found from this file's own location; a section run from the
% editor has no mfilename, so fall back to the file open in the editor.
fabric_dir = fileparts(mfilename('fullpath'));
if isempty(fabric_dir) || ~isfile(fullfile(fabric_dir,'fabric_task.m'))
  try
    fabric_dir = fileparts(matlab.desktop.editor.getActiveFilename);
  catch
    fabric_dir = '';
  end
end
if isempty(fabric_dir) || ~isfile(fullfile(fabric_dir,'fabric_task.m'))
  error('run_fabric:paths', ['Cannot locate the opr_fabric directory. ' ...
    'Run this file as run(''<clone>/opr_fabric/run_fabric.m'').']);
end
addpath(fileparts(fabric_dir), fabric_dir);
% test/stubs shadows OPR's opr_* functions for the offline tests and the
% scratch batches; if one of those ran in this session, take it back off
if contains([pathsep path pathsep], [pathsep fullfile(fabric_dir,'test','stubs') pathsep])
  rmpath(fullfile(fabric_dir,'test','stubs'));
end
global gRadar;
if isempty(gRadar) || ~exist('read_param_xls','file')
  error('run_fabric:opr', ['OPR is not set up in this MATLAB session ' ...
    '(gRadar is empty). Run your OPR startup.m first.']);
end

%% User Setup
% =====================================================================
param_override = [];

if 1
  params = read_param_xls(opr_filename_param('accum_param_2024_Antarctica_Ground2.xlsx'));

  % Example to run specific segments and frames by overriding parameter spreadsheet values
  params = opr_set_params(params,'cmd.generic',0);
  params = opr_set_params(params,'cmd.generic',1,'day_seg','20250108_02');
  params = opr_set_params(params,'cmd.frms',[9]);

  % Input product: must match polarimetric.out_path used in run_polarimetric.
  % For this season the SNAPHU-unwrapped product is Lilien's
  % CSARP_polarimetric_unwrap_dlilien (frames 1, 3, 7, 9 of 20250108_02)
  params = opr_set_params(params,'fabric.in_path','polarimetric_unwrap_dlilien');
  % Output product: fabric_joint keeps the CSARP_fabric stripping outputs
  % untouched for side-by-side comparison
  params = opr_set_params(params,'fabric.out_path','fabric_joint');

  param_override.fabric.frm_types = {0,0,-1,-1,-1}; % Only do frames that are SAR processed

  % Traveltime-difference estimation
  params = opr_set_params(params,'fabric.coherence_threshold',0.5);
  params = opr_set_params(params,'fabric.use_snaphu_phase',true);
  params = opr_set_params(params,'fabric.blend_coreg_en',true);
  params = opr_set_params(params,'fabric.phase_sign',0); % 0 = auto from coregistration
  params = opr_set_params(params,'fabric.block_size',1000);

  % Inversion setup
  params = opr_set_params(params,'fabric.inversion','joint'); % smoothness-regularized joint solve ('stripping' for the exact per-interval solve)
  params = opr_set_params(params,'fabric.num_intervals',10);
  params = opr_set_params(params,'fabric.half_offset',0); % half Tx-Rx separation [m]
  % Firn-ice column model (see ptt.defaultParams). H = local ice thickness;
  % bco_depth = bubble close-off depth below surface [m]; lam_z profile is
  % ASSUMED (isotropic-vertical default is fine near the surface).
  params = opr_set_params(params,'fabric.ptt',struct( ...
    'H',2000,'bco_depth',60,'lam_z_sfc',1/3,'lam_z_bed',1/3));

end

if ~batchStartupOptionUsed
  dbstop if error; % not under matlab -batch: it would park the run at K>> forever
end
if 0
  param_override.cluster.type = 'slurm';
  param_override.cluster.max_jobs_active = 128;
elseif 1
  param_override.cluster.type = 'debug';
end
% param_override.cluster.rerun_only = true;
% run_chain: run the chain here. false = OPR's usual save-only behavior;
% then run it later with cluster_run(cluster_load_chain(<chain_id>))
run_chain = true;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~opr_generic_en(param)
    continue;
  end
  ctrl_chain{end+1} = fabric(param,param_override);
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

if run_chain
  ctrl_chain = cluster_run(ctrl_chain);
end
