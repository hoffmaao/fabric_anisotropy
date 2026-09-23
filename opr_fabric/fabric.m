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
if ~ischar(param.fabric.dtau_source) ...
    || ~any(strcmp(param.fabric.dtau_source,{'phase','coreg','deltak'}))
  error('param.fabric.dtau_source must be ''phase'', ''coreg'' or ''deltak''.');
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

% seam_mask_en: exclude the fast-time samples straddling the waveform
% image-combination boundaries, where the trace is a weighted blend of a
% short and a long pulse (different bandwidth, gain and system delay) and
% the interferometric phase is not an ice property. The boundaries come
% from the product's own param.array.img_comb, so they follow the frame
% (0.9 us / 0.1 us blend on the accum3 grids, 8 us / 1 us on the deeper
% settings, none on single-image frames). See ptt.imgCombSeam.
if ~isfield(param.fabric,'seam_mask_en') || isempty(param.fabric.seam_mask_en)
  param.fabric.seam_mask_en = true;
end

% seam_mask_win: guard width as a multiple of the blend window each
% boundary declares. The masked band is asymmetric - OPR crossfades over
% [t_comb, t_comb+window], so the contamination runs forward - and spans
% [t_comb - k*window, t_comb + (1+k)*window] for k = seam_mask_win.
if ~isfield(param.fabric,'seam_mask_win') || isempty(param.fabric.seam_mask_win)
  param.fabric.seam_mask_win = 1;
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

% ref_band_twtt: width of the band below ref_twtt_offset that the
% reference is AVERAGED over (ptt.blendTraveltime), and that the inversion
% nodes must sit below (ptt.invertBlocks). A single-bin reference injects
% its own error as a constant into every node, and a constant can only
% land in the shallowest interval's dlam - the mechanism behind the
% spurious non-zero near-surface fabric. The residual after band
% averaging is carried by the joint inversion's reference-offset
% nuisance, and the first interval is flagged ref_degenerate.
if ~isfield(param.fabric,'ref_band_twtt') || isempty(param.fabric.ref_band_twtt)
  param.fabric.ref_band_twtt = 100e-9;
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
% mcc runs in a fresh process that sees only its own startup.m path, not
% this session's, so a checkout that run_fabric.m put on the path is
% invisible to it: the compiled job then dies on the cluster with
% "Unable to resolve the name 'ptt.blendTraveltime'". Name this checkout
% explicitly - the task by full path, -I for its folder and the repo root,
% and -a for the whole +ptt package. The flags ride in hidden_depend_funs
% at check level 0, which cluster_compile appends to the mcc command as-is.
% A fabric_task entry left in hidden_depend_funs by the old deployment
% steps would make mcc resolve a second, possibly stale, copy by name, so
% drop it: only this checkout's task goes in, by full path.
fabric_dir = fileparts(mfilename('fullpath'));
repo_dir = fileparts(fabric_dir);
mcc_flags = {sprintf('-I ''%s'' -I ''%s'' -a ''%s''', fabric_dir, repo_dir, ...
  fullfile(repo_dir,'+ptt')), 0};
hidden_depend_funs = ctrl.cluster.hidden_depend_funs;
if ~isempty(hidden_depend_funs)
  keep = true(size(hidden_depend_funs));
  for dep_idx = 1:numel(hidden_depend_funs)
    [~,dep_name] = fileparts(hidden_depend_funs{dep_idx}{1});
    keep(dep_idx) = ~strcmp(dep_name,'fabric_task');
  end
  hidden_depend_funs = hidden_depend_funs(keep);
end
% Always compile under slurm/torque: cluster_job is shared by every OPR
% task, so another task's recompile drops fabric_task and +ptt from it, and
% the date check cannot see +ptt functions reached through handles.
force_compile = ctrl.cluster.force_compile ...
  || any(strcmpi(ctrl.cluster.type,{'slurm','torque'}));
cluster_compile({fullfile(fabric_dir,'fabric_task.m')}, ...
  [hidden_depend_funs {mcc_flags}],force_compile,ctrl);

ctrl_chain = {};

% Input product directory, resolved exactly as fabric_task.m resolves it.
% A frame whose product is missing would otherwise run a task that warns,
% writes nothing and still completes - so check here, before submitting.
in_dir = opr_filename_out(param,param.fabric.in_path);
num_in = 0; num_missing = 0;

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

  if param.fabric.img == 0
    in_fn = fullfile(in_dir, sprintf('Data_%s.mat', frm_id));
  else
    in_fn = fullfile(in_dir, sprintf('Data_img_%02d_%s.mat', param.fabric.img, frm_id));
  end
  if ~exist(in_fn,'file')
    fprintf('Skipping %s (no input product %s)\n', frm_id, in_fn);
    num_missing = num_missing + 1;
    continue;
  end
  num_in = num_in + 1;

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
  if strcmp(param.fabric.dtau_source,'deltak')
    % The ladder adds two full-frame range FFTs, an ifft pair per sub-band
    % (12+4+2 bands) and three conv2 passes on top of the normal chain
    dparam.cpu_time = 7200;
    dparam.mem = 24e9;      % holds ref/sec SLC spectra + band products
  else
    dparam.cpu_time = 3600; % Much lighter than polarimetric_task
    dparam.mem = 12e9;
  end

  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);

end

if num_in == 0 && num_missing > 0
  % Nothing to do because fabric.in_path is wrong, not because the work is
  % done: say which products this season does have. An absolute in_path is
  % outside the season tree, so there is no season to list.
  in_path = param.fabric.in_path;
  if in_path(1) == '/' || in_path(1) == '\' ...
      || (ispc && (contains(in_path,':\') || contains(in_path,':/')))
    avail_msg = '';
  else
    season_dir = fileparts(fileparts(in_dir));
    avail = dir(fullfile(season_dir,'CSARP_polarimetric*'));
    if isempty(avail)
      avail_msg = sprintf(' No CSARP_polarimetric* products in\n  %s', ...
        season_dir);
    else
      avail_msg = sprintf(' Polarimetric products in this season:%s', ...
        sprintf('\n  %s', avail.name));
    end
  end
  error('fabric:noInput', ['No requested frame of %s has an input product ' ...
    'in\n  %s\nCheck fabric.in_path.%s'], param.day_seg, in_dir, avail_msg);
end

ctrl = cluster_save_dparam(ctrl);

ctrl_chain{end+1} = ctrl;

fprintf('Done %s\n', datestr(now));
