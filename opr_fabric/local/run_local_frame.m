%RUN_LOCAL_FRAME Single-frame polarimetric + fabric test on local data.
%   Runs the real OPR polarimetric_task and our fabric_task directly (no
%   cluster, no frames/records support files) on one quad-pol frame
%   downloaded from the public portal with scripts/fetch_frame.sh.
%
%   Run this in MATLAB (the coregistration code uses an arguments block
%   that Octave cannot parse), e.g. the browser MATLAB of the
%   fabric-matlab container: cd to this directory and run run_local_frame.
%
%   The gRadar-dependent OPR path helpers are shadowed by the stubs in
%   ../test/stubs, which map CSARP_* products under data_root. Outputs
%   land in <data_root>/CSARP_polarimetric and <data_root>/CSARP_fabric.

%% Paths (container defaults; falls back to host-side locations)
opr_root = '/home/matlab/opr/matlab';
data_root = '/home/matlab/data/accum/2024_Antarctica_Ground2';
if ~exist(opr_root,'dir')
  opr_root = fullfile(getenv('HOME'),'projects','opr','matlab');
  data_root = fullfile(getenv('HOME'),'data','opr','accum','2024_Antarctica_Ground2');
end

this_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));

addpath(genpath(opr_root));                       % OPR toolbox
addpath(proj_root);                               % +ptt
addpath(fullfile(proj_root,'opr_fabric'));        % fabric_task
addpath(fullfile(proj_root,'opr_fabric','test','stubs')); % shadow opr_* path helpers

%% Frame selection (matches scripts/fetch_frame.sh defaults)
param = [];
param.day_seg = '20250108_02';
param.season_name = '2024_Antarctica_Ground2';
param.radar_name = 'accum3';
param.load.frm = 9;
param.opr_file_lock = false;
param.stub_out_root = data_root;

%% Stage 1: polarimetric_task (real OPR code, all defaults set explicitly)
pp = [];
pp.HH_path = 'standard_HH';
pp.VV_path = 'standard_VV';
pp.HV_path = 'standard_HV';
pp.VH_path = 'standard_VH';
pp.img = 0;
pp.out_path = 'polarimetric';
pp.out_file_exts = {'.jpg'};
pp.synth_rot_deg = 0;            % rotate basis toward principal axes later
pp.mlook_window = [5 15];
pp.power_mlook_window = [1 15];
pp.min_rbin = 100;               % values John used for this segment
pp.max_rbin = 6700;              % (contiguous internal layers)
pp.chan_equal_HH = 1; pp.chan_equal_VV = 1;
pp.chan_equal_HV = 1; pp.chan_equal_VH = 1;
pp.rotation_movie.en = false;    % skip the (slow) rotation movie
pp.coregistration.en = true;
pp.coregistration.use_results_en = true;
pp.coregistration.Tt = 51;
pp.coregistration.Tx = 101;
pp.coregistration.overlap_t = 25;
pp.coregistration.overlap_x = 25;
pp.coregistration.search_t = 5;
pp.coregistration.search_x = 5;
pp.coregistration.one_dim_search_en = true;
pp.snaphu_en = false;            % no snaphu binary locally; the fabric
                                 % blend then unwraps via coregistration
pp.snaphu_path = '';
param.polarimetric = pp;

fprintf('=== Stage 1: polarimetric_task ===\n');
success = polarimetric_task(param);
assert(success,'polarimetric_task failed');

%% Stage 2: fabric_task
pf = [];
pf.in_path = 'polarimetric';
pf.out_path = 'fabric';
pf.img = 0;
pf.out_file_exts = {'.jpg'};
pf.fc = 750e6;                   % accum3 600-900 MHz; refine from wfs
pf.dtau_source = 'coreg';        % public CSARP_standard_* files are
                                 % detected POWER (no phase): coregistration
                                 % offsets are the only valid dtau source.
                                 % Switch to 'phase' when running on
                                 % standardphase_* products from the servers
pf.phase_sign = 0;               % auto from coregistration ('phase' mode)
pf.use_snaphu_phase = false;     % wrapped-phase + per-pixel fringe blend
pf.blend_coreg_en = true;
pf.coherence_threshold = 0.5;
pf.min_coverage = 0.3;
pf.block_size = 1000;
pf.inversion = 'joint';          % same solver as run_fabric_on_product
pf.num_intervals = 8;
pf.ref_twtt_offset = 50e-9;
pf.half_offset = 0;              % co-located crossed bowties (switch-gated):
                                 % phase centers coincide, true nadir geometry
pf.ptt = struct('H', 2000, ...   % standing default; H cancels in the
  'bco_depth', 60, ...           % twtt-depth map at nadir (only bco_depth
  'lam_z_sfc', 1/3, 'lam_z_bed', 1/3); % matters) and only gates the base cutoff
param.fabric = pf;

fprintf('=== Stage 2: fabric_task ===\n');
success = fabric_task(param);
assert(success,'fabric_task failed');

fprintf('\nOutputs:\n  %s\n  %s\n', ...
  fullfile(data_root,'CSARP_polarimetric',param.day_seg), ...
  fullfile(data_root,'CSARP_fabric',param.day_seg));
