function ctrl_chain = fabric(param,param_override)
% ctrl_chain = fabric(param,param_override)
%
% Infers the horizontal ice-fabric contrast dlam = lam_x - lam_y versus
% depth and along-track position from the CSARP_polarimetric product
% (polarimetric.m/polarimetric_task.m output). The traveltime difference
% between the two synthesized polarizations is estimated by blending the
% coregistration row offsets (unambiguous, coarse) with the multilooked
% interferogram phase (precise, 2*pi-ambiguous; SNAPHU-unwrapped when
% available), then inverted per along-track block (smoothness-regularized
% joint solve by default, or exact layer stripping; param.fabric.inversion)
% through the Maxwell-Garnett firn model of Rathmann (2026),
% doi:10.1098/rspa.<pending> (see the +ptt package, which must be on the
% MATLAB path).
%
% An output .mat file (CSARP_fabric) and image files are output.
%
% See also: run_fabric.m, fabric_task.m, polarimetric.m, +ptt

%% General Setup
% =====================================================================
if exist('param_override','var')
  param = merge_structs(param, param_override);
end
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input arguments check and setup
% =========================================================================

% er_ice, c = speed of light
physical_constants;
[output_dir,radar_type,radar_name] = opr_output_dir(param.radar_name);

% cmd worksheet inputs
% -------------------------------------------------------------------------

% Load the frames file
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

% fabric worksheet inputs
% -------------------------------------------------------------------------

% in_path: CSARP_polarimetric product to read (polarimetric.out_path)
if ~isfield(param.fabric,'in_path') || isempty(param.fabric.in_path)
  param.fabric.in_path = 'polarimetric';
end

if ~isfield(param.fabric,'img') || isempty(param.fabric.img)
  param.fabric.img = 0;
end

if ~isfield(param.fabric,'out_path') || isempty(param.fabric.out_path)
  param.fabric.out_path = 'fabric';
end

if ~isfield(param.fabric,'out_file_exts') || isempty(param.fabric.out_file_exts)
  param.fabric.out_file_exts = {'.jpg'};
end

% fc: center frequency [Hz] used to convert interferogram phase to
% traveltime. Default: derived from param.radar.wfs(1).
if ~isfield(param.fabric,'fc') || isempty(param.fabric.fc)
  param.fabric.fc = (param.radar.wfs(1).f0 + param.radar.wfs(1).f1)/2;
end

% phase_sign: sign s in dtau = s*phase/(2*pi*fc). 0 (default) estimates the
% sign per frame by regressing the phase against the coregistration row
% offsets; +/-1 forces it.
if ~isfield(param.fabric,'phase_sign') || isempty(param.fabric.phase_sign)
  param.fabric.phase_sign = 0;
end

% use_snaphu_phase: use snaphu_out_phase from the polarimetric product when
% present (recommended); otherwise the wrapped interferogram phase is used
% and fringe ambiguities are resolved per pixel from the coregistration.
if ~isfield(param.fabric,'use_snaphu_phase') || isempty(param.fabric.use_snaphu_phase)
  param.fabric.use_snaphu_phase = true;
end

% dtau_source: 'phase' (default; needs phase-preserved input products like
% standardphase_*), 'coreg' (traveltime differences from coregistration
% row offsets alone; the only valid choice when the polarimetric product
% was formed from detected-power echograms such as the public
% CSARP_standard_* files, where interferogram phase is meaningless), or
% 'deltak' (split-spectrum ladder over the ref/sec SLC spectra: absolute
% dtau per pixel with NO phase unwrapping and no fringe blending; needs
% the polarimetric product saved with ref and the unregistered sec, and
% more task memory; see ptt.deltakTraveltime and param.fabric.deltak)
if ~isfield(param.fabric,'dtau_source') || isempty(param.fabric.dtau_source)
  param.fabric.dtau_source = 'phase';
end

% blend_coreg_en: resolve integer-fringe offsets of the phase-derived
% traveltime differences using the coregistration row offsets
if ~isfield(param.fabric,'blend_coreg_en') || isempty(param.fabric.blend_coreg_en)
  param.fabric.blend_coreg_en = true;
end

% coherence_threshold: pixels with multilooked coherence below this are
% excluded from block averages
if ~isfield(param.fabric,'coherence_threshold') || isempty(param.fabric.coherence_threshold)
  param.fabric.coherence_threshold = 0.5;
end

% min_coverage: minimum fraction of coherent pixels for a (bin, block) to
% contribute to the averaged dtau profile
if ~isfield(param.fabric,'min_coverage') || isempty(param.fabric.min_coverage)
  param.fabric.min_coverage = 0.3;
end

% block_size: number of range lines averaged into one inversion
if ~isfield(param.fabric,'block_size') || isempty(param.fabric.block_size)
  param.fabric.block_size = 1000;
end

% num_intervals: number of depth intervals (piecewise-constant dlam nodes)
if ~isfield(param.fabric,'num_intervals') || isempty(param.fabric.num_intervals)
  param.fabric.num_intervals = 10;
end

% inversion: 'stripping' (exact per-interval layer stripping) or 'joint'
% (smoothness-regularized joint solve; robust to noisy dtau increments
% that make stripping oscillate and peg the eigenvalue bounds)
if ~isfield(param.fabric,'inversion') || isempty(param.fabric.inversion)
  param.fabric.inversion = 'joint';
end

% reg: regularization strength for the joint inversion (dimensionless
% relative to the mean data sensitivity; see invertHorizontalFabricJoint)
if ~isfield(param.fabric,'reg') || isempty(param.fabric.reg)
  param.fabric.reg = 0.05;
end

% ref_twtt_offset: two-way time below the surface return where dtau is
% referenced to zero (avoids surface sidelobes; no birefringence above)
if ~isfield(param.fabric,'ref_twtt_offset') || isempty(param.fabric.ref_twtt_offset)
  param.fabric.ref_twtt_offset = 50e-9;
end

% half_offset: half the Tx-Rx antenna separation [m]. 0 for the ground
% accum radar: switch-gated crossed bowties share one phase center.
if ~isfield(param.fabric,'half_offset') || isempty(param.fabric.half_offset)
  param.fabric.half_offset = 0;
end

% ptt: firn-ice column model parameters (see ptt.defaultParams). H is the
% ice-column thickness; at nadir the twtt-depth mapping depends only on
% bco_depth (H cancels), so H just needs to exceed the coherent depth
% range and 2000 m is the standing default. bco_depth [m below surface]
% overrides zhat_bco. lam_z_sfc/lam_z_bed set the ASSUMED vertical
% eigenvalue profile (not constrained by common-offset data).
if ~isfield(param.fabric,'ptt') || isempty(param.fabric.ptt)
  param.fabric.ptt = struct();
end
if ~isfield(param.fabric.ptt,'H') || isempty(param.fabric.ptt.H)
  param.fabric.ptt.H = 2000;
end
if ~isfield(param.fabric.ptt,'bco_depth') || isempty(param.fabric.ptt.bco_depth)
  param.fabric.ptt.bco_depth = 60;
end

if ~isfield(param.fabric,'frm_types') || isempty(param.fabric.frm_types)
  param.fabric.frm_types = {0,0,-1,-1,-1}; % Default to only doing frames with SAR processing enabled and that are not disabled in general
end

%% Setup cluster
% =====================================================================
ctrl = cluster_new_batch(param);
cluster_compile({'fabric_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);

ctrl_chain = {};

%% Block: Create tasks (one for each frame)
% =====================================================================
sparam.argsin{1} = param; % Static parameters
sparam.task_function = 'fabric_task';
sparam.num_args_out = 1;

for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);

  % Check proc_mode from frames file that contains this frames type and
  % make sure the user has specified to process this frame type
  if opr_proc_frame(frames.proc_mode(frm),param.fabric.frm_types)
    fprintf('%s %s_%03i (%i of %i) (%s)\n', sparam.task_function, param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  frm_id = sprintf('%s_%03d', param.day_seg, frm);

  % Prepare task inputs
  % =================================================================
  dparam = [];
  dparam.argsin{1}.load.frm = frm;

  % Create success condition
  % =================================================================
  dparam.file_success = {};

  % .mat file
  out_dir = opr_filename_out(param,param.fabric.out_path);
  if param.fabric.img == 0
    fn_name = sprintf('Data_%s.mat', frm_id);
  else
    fn_name = sprintf('Data_img_%02d_%s.mat', param.fabric.img, frm_id);
  end

  dparam.file_success{end+1} = fullfile(out_dir,fn_name);
  if ~ctrl.cluster.rerun_only && exist(dparam.file_success{end},'file')
    delete(dparam.file_success{end});
  end

  % image files (e.g. {'.jpg', '.fig'})
  [~,fn_name,fn_ext] = fileparts(fn_name);
  for file_ext = param.fabric.out_file_exts
    file_ext = file_ext{1};

    dparam.file_success{end+1} = fullfile(out_dir,sprintf('%s_dlam%s',fn_name,file_ext));
    if ~ctrl.cluster.rerun_only && exist(dparam.file_success{end},'file')
      delete(dparam.file_success{end});
    end
    dparam.file_success{end+1} = fullfile(out_dir,sprintf('%s_dtau%s',fn_name,file_ext));
    if ~ctrl.cluster.rerun_only && exist(dparam.file_success{end},'file')
      delete(dparam.file_success{end});
    end
  end

  % Rerun only mode: Test to see if we need to run this task
  % =================================================================
  dparam.notes = sprintf('%s:%s:%s:%s %s_%03d (%d of %d)', ...
    sparam.task_function, param.radar_name, param.season_name, out_dir, param.day_seg, frm, frm_idx, length(param.cmd.frms));
  if ctrl.cluster.rerun_only
    if ~cluster_file_success(dparam.file_success)
      fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
        dparam.notes, datestr(now));
      continue;
    end
  end

  % Create task
  % =================================================================
  dparam.cpu_time = 3600; % Much lighter than polarimetric_task
  if strcmp(param.fabric.dtau_source,'deltak')
    dparam.mem = 24e9;    % holds ref/sec SLC spectra + band products
  else
    dparam.mem = 12e9;
  end

  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);

end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain{end+1} = ctrl;

fprintf('Done %s\n', datestr(now));
