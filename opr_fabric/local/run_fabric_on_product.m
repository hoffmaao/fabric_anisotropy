%RUN_FABRIC_ON_PRODUCT Fabric inversion on a real CSARP_polarimetric file.
%   Runs fabric_task on a downloaded polarimetric product (e.g. the
%   extracted CSARP_polarimetric_unwrap_dlilien frame from mem1). Unlike
%   run_local_frame.m this needs no MATLAB-only toolbox code, so it runs
%   in Octave or MATLAB.
%
%   Expects the product at
%   <data_root>/CSARP_polarimetric_unwrap/<day_seg>/Data_<day_seg>_FFF.mat

%% Paths
data_root = '/data/accum/2024_Antarctica_Ground2'; % Octave container mount
if ~exist(data_root,'dir')
  data_root = fullfile(getenv('HOME'),'data','opr','accum','2024_Antarctica_Ground2');
end

this_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));
addpath(proj_root);                               % +ptt
addpath(fullfile(proj_root,'opr_fabric'));        % fabric_task
addpath(fullfile(proj_root,'opr_fabric','test','stubs')); % opr_* stubs

%% Param setup
param = [];
param.day_seg = '20250108_02';
param.season_name = '2024_Antarctica_Ground2';
param.radar_name = 'accum3';
param.load.frm = 9;
param.opr_file_lock = false;
param.stub_out_root = data_root;

pf = [];
pf.in_path = 'polarimetric_unwrap';
pf.out_path = 'fabric_joint';
pf.img = 0;
pf.out_file_exts = {'.png'};
pf.fc = 750e6;                  % deep waveform 600-900 MHz
pf.dtau_source = 'phase';       % real standardphase-based interferogram
pf.use_snaphu_phase = true;     % SNAPHU output present in this product
pf.blend_coreg_en = true;
pf.phase_sign = 0;              % auto from coregistration
pf.coherence_threshold = 0.5;
pf.min_coverage = 0.3;
pf.block_size = 1000;           % 4346 rlines -> 5 blocks
pf.num_intervals = 10;
pf.ref_twtt_offset = 50e-9;
pf.half_offset = 0;             % colocated crossed bowties
pf.inversion = 'joint';
pf.reg = 0.05;
pf.ptt = struct('H', 2000, 'bco_depth', 60, 'lam_z_sfc', 1/3, 'lam_z_bed', 1/3);
param.fabric = pf;

%% Run
success = fabric_task(param);
assert(success, 'fabric_task failed');

out_fn = fullfile(data_root,['CSARP_' pf.out_path],param.day_seg, ...
  sprintf('Data_%s_%03d.mat',param.day_seg,param.load.frm));
out = load(out_fn);
fprintf('\nBlocks: %d, intervals: %d\n', size(out.dlam,2), size(out.dlam,1));
fprintf('Phase sign: %+d, fringes: %s\n', out.phase_sign, mat2str(out.blend_fringes));
fprintf('\nBlock-median dlam profile (depth is interval bottom, m):\n');
for k = 1:size(out.dlam,1)
  fprintf('  %6.1f m   dlam = %+7.4f   (coh %.2f)\n', ...
    median(out.dlam_bot_depth(k,:),'omitnan'), ...
    median(out.dlam(k,:),'omitnan'), median(out.dlam_quality(k,:),'omitnan'));
end
