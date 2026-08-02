% script run_fabric
%
% Script for running fabric.m (usually just used for debugging).
%
% Requires: the +ptt package on the MATLAB path, and CSARP_polarimetric
% products from polarimetric.m (run with coregistration.en = true and
% ideally snaphu_en = true).
%
% See also: run_master.m, master.m, fabric.m, fabric_task.m,
% run_polarimetric.m, polarimetric.m

%% User Setup
% =====================================================================
param_override = [];

if 1
  params = read_param_xls(opr_filename_param('accum_param_2024_Antarctica_Ground2.xlsx'));

  % Example to run specific segments and frames by overriding parameter spreadsheet values
  params = opr_set_params(params,'cmd.generic',0);
  params = opr_set_params(params,'cmd.generic',1,'day_seg','20250108_02');
  params = opr_set_params(params,'cmd.frms',[9]);

  % Input product: must match polarimetric.out_path used in run_polarimetric
  params = opr_set_params(params,'fabric.in_path','polarimetric_unwrap');
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

dbstop if error;
if 0
  param_override.cluster.type = 'slurm';
  param_override.cluster.max_jobs_active = 128;
elseif 1
  param_override.cluster.type = 'debug';
end
% param_override.cluster.rerun_only = true;

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
