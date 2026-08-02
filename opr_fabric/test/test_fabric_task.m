%TEST_FABRIC_TASK End-to-end test of fabric_task on a synthetic product.
%   Builds a synthetic CSARP_polarimetric frame from a known fabric using
%   the ptt forward model: interferogram phase (with matched-filter sign,
%   an arbitrary unwrapping constant, channel timing bias, and noise),
%   coherence that decays with depth, and noisy coregistration row offsets.
%   Then runs the real fabric_task (with OPR support functions stubbed) and
%   compares the inferred dlam profile against the truth.

clear;
rng(3);
t0 = tic;

thisDir = fileparts(mfilename('fullpath'));
projRoot = fullfile(thisDir, '..', '..');
addpath(projRoot);                    % +ptt
addpath(fullfile(thisDir, '..'));     % fabric_task.m
addpath(fullfile(thisDir, 'stubs'));  % OPR support stubs

outRoot = fullfile(thisDir, 'output');
if exist(outRoot, 'dir'), rmdir(outRoot, 's'); end

%% True column: girdle strengthening over the top ~300 m
parT = ptt.defaultParams();
parT.H = 2000;
parT.zhat_bco = 1 - 60/parT.H; % bubble close-off at 60 m depth
parT.lam_x_sfc = 1/3; parT.lam_z_sfc = 1/3;
% Linear-in-zhat endpoints chosen so dlam ~ -0.3 at 300 m depth
parT.lam_x_bed = 1/3 - 0.5*2000/300*0.3; parT.lam_z_bed = 1/3;
% lam_x_bed < -1 is kept as configured; only the top 300 m is observed

%% Synthetic polarimetric product
fc = 750e6;
dt = 2e-9;
Nt = 2200; Nx = 3000;
Time = (0:Nt-1).'*dt;             % 0 to 4.4 us
surf_twtt = 0.3e-6;
Surface = surf_twtt*ones(1,Nx);

% Vertical dtau(t) profile from the model (t_x - t_y, seconds)
zhat_fine = linspace(1, 1-400/parT.H, 4001).';
P = ptt.columnProfiles(parT, zhat_fine);
depth_fine = parT.H*(1-zhat_fine);
twtt_fine = 2e-9*cumtrapz(depth_fine, P.S(:,2));
dtau_fine = 2e-9*cumtrapz(depth_fine, P.S(:,1) - P.S(:,2)); % t_x - t_y [s]
t_rel = Time - surf_twtt;
dtau_true = interp1(twtt_fine, dtau_fine, t_rel, 'linear', 0);
dtau_true(t_rel < 0) = 0;
dtau_map = repmat(dtau_true, 1, Nx);

% Coherence: high below surface, decaying to nothing below ~300 m depth
depth_map = interp1(twtt_fine, depth_fine, t_rel, 'linear', 'extrap');
coh_prof = 0.9 ./ (1 + exp((depth_map - 300)/25));
coh_prof(t_rel < 0) = 0.05;   % no coherence above the surface
coherence = max(0.01, repmat(coh_prof, 1, Nx) + 0.03*randn(Nt, Nx));

% Interferogram phase: matched-filter convention phi = -2*pi*fc*dtau, plus
% a channel timing bias (removed by surface referencing), phase noise where
% coherence is low, and an arbitrary unwrapping constant in snaphu output
timing_bias = 0.8e-9;
phase_noise = 0.4*randn(Nt, Nx) ./ max(coherence, 0.05) * 0.1;
phase_true = -2*pi*fc*(dtau_map + timing_bias);
snaphu_out_phase = phase_true + phase_noise + 2*pi*3; % unwrap constant
interferogram_mlook = single(coherence .* exp(1i*(phase_true + phase_noise)));
interferogram_coherence = single(coherence);

% Coregistration row offsets: dtau/dt plus coarse noise (0.3 bins)
row_offset = (dtau_map + timing_bias)/dt + 0.3*randn(Nt, Nx);

GPS_time = 1e9 + (0:Nx-1)*0.1;
Latitude = -75 + (0:Nx-1)*1e-5;
Longitude = 123 + (0:Nx-1)*1e-5;
Elevation = 3000*ones(1,Nx);
Bottom = nan(1,Nx);
param_records = struct('note','synthetic');
param_polarimetric = struct('note','synthetic');
file_type = 'polarimetric';
file_version = '1';

day_seg = '20250108_02';
in_dir = fullfile(outRoot, 'CSARP_polarimetric_unwrap', day_seg);
mkdir(in_dir);
save(fullfile(in_dir, sprintf('Data_%s_009.mat', day_seg)), ...
  'interferogram_mlook','interferogram_coherence','snaphu_out_phase', ...
  'row_offset','Time','GPS_time','Latitude','Longitude','Elevation', ...
  'Surface','Bottom','param_records','param_polarimetric', ...
  'file_type','file_version');

%% Assemble the task param struct (what fabric.m + cluster would provide)
param = [];
param.day_seg = day_seg;
param.radar_name = 'accum3';
param.season_name = '2024_Antarctica_Ground2';
param.load.frm = 9;
param.opr_file_lock = false;
param.stub_out_root = outRoot;
param.radar.wfs(1).f0 = 600e6;
param.radar.wfs(1).f1 = 900e6;

param.fabric.in_path = 'polarimetric_unwrap';
param.fabric.out_path = 'fabric';
param.fabric.img = 0;
param.fabric.out_file_exts = {'.png'};
param.fabric.fc = fc;
param.fabric.phase_sign = 0;         % exercise auto sign detection
param.fabric.use_snaphu_phase = true;
param.fabric.blend_coreg_en = true;
param.fabric.coherence_threshold = 0.5;
param.fabric.min_coverage = 0.3;
param.fabric.block_size = 1000;
param.fabric.num_intervals = 8;
param.fabric.ref_twtt_offset = 50e-9;
param.fabric.half_offset = 0;
param.fabric.ptt = struct('H', parT.H, 'bco_depth', 60, ...
  'lam_z_sfc', 1/3, 'lam_z_bed', 1/3);

%% Run the task
success = fabric_task(param);
assert(success, 'fabric_task did not succeed');

%% Compare against truth
out = load(fullfile(outRoot, 'CSARP_fabric', day_seg, ...
  sprintf('Data_%s_009.mat', day_seg)));

fprintf('\nBlock 1: interval (depth m)   true dlam   inferred\n');
max_err = 0;
for k = 1:size(out.dlam,1)
  dmid = (out.dlam_top_depth(k,1) + out.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  dlam_true_k = lam_mid.lam(1) - lam_mid.lam(2);
  err = abs(out.dlam(k,1) - dlam_true_k);
  max_err = max(max_err, err);
  fprintf('  %5.0f - %5.0f          %8.3f  %8.3f\n', ...
    out.dlam_top_depth(k,1), out.dlam_bot_depth(k,1), dlam_true_k, out.dlam(k,1));
end
fprintf('Max |error| block 1: %.3f\n', max_err);
fprintf('Phase sign detected: %+d (expected -1)\n', out.phase_sign);
fprintf('Blend fringes per block: %s (expected 0; constants removed by referencing)\n', mat2str(out.blend_fringes));

assert(out.phase_sign == -1, 'auto sign detection failed');
assert(max_err < 0.05, 'inferred dlam deviates from truth by %.3f', max_err);
nblk_ok = sum(all(isfinite(out.dlam),1));
fprintf('Blocks fully inverted: %d of %d\n', nblk_ok, size(out.dlam,2));
assert(nblk_ok == size(out.dlam,2), 'not all blocks inverted');

%% Coregistration-only mode (detected-power products: no usable phase)
param.fabric.dtau_source = 'coreg';
param.fabric.out_path = 'fabric_coreg';
success = fabric_task(param);
assert(success, 'fabric_task (coreg mode) did not succeed');

outc = load(fullfile(outRoot, 'CSARP_fabric_coreg', day_seg, ...
  sprintf('Data_%s_009.mat', day_seg)));
max_err_c = 0;
for k = 1:size(outc.dlam,1)
  dmid = (outc.dlam_top_depth(k,1) + outc.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  max_err_c = max(max_err_c, abs(outc.dlam(k,1) - (lam_mid.lam(1) - lam_mid.lam(2))));
end
fprintf('Coreg-only mode: max |error| block 1: %.3f\n', max_err_c);
assert(max_err_c < 0.1, 'coreg-only dlam deviates from truth by %.3f', max_err_c);

fprintf('\nPASS (%.1f s)\n', toc(t0));
