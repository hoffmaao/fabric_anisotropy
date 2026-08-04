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

% Waveform image-combination seam: OPR crossfades two waveform images over
% [t_comb, t_comb+window], where the trace is a blend of a short and a long
% pulse rather than an ice property. Declare it the way a real product does
% (param.array.img_comb = [t_comb mult window], mult = -Inf pinning the
% boundary to a fixed traveltime) and contaminate the crossfade to match: a
% bogus group delay and a coherence step that still clears the coherence
% threshold, so only ptt.imgCombSeam can remove it.
seam_t_comb = 1.5e-6;
seam_win = 0.1e-6;
seam_bias = 3e-9;                 % [s] spurious delay across the blend
in_seam = Time >= seam_t_comb & Time <= seam_t_comb + seam_win;
dtau_map(in_seam,:) = dtau_map(in_seam,:) + seam_bias;
coherence(in_seam,:) = 0.6;

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
param_polarimetric = struct('note','synthetic', ...
  'array', struct('img_comb', [seam_t_comb -Inf seam_win]));
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

% Intervals bounded by a node whose dtau was interpolated across the seam
% gap are fabricated, not measured (that is what dlam_interpolated marks),
% and layer stripping carries the fabricated increment into the interval
% below it as well; both are excluded from the accuracy check.
seam_gap = @(f,k) f(k,1) == 1 || (k > 1 && f(k-1,1) == 1);

fprintf('\nBlock 1: interval (depth m)   true dlam   inferred\n');
max_err = 0;
for k = 1:size(out.dlam,1)
  dmid = (out.dlam_top_depth(k,1) + out.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  dlam_true_k = lam_mid.lam(1) - lam_mid.lam(2);
  if seam_gap(out.dlam_interpolated, k)
    fprintf('  %5.0f - %5.0f          %8.3f  %8.3f  (seam gap)\n', ...
      out.dlam_top_depth(k,1), out.dlam_bot_depth(k,1), dlam_true_k, out.dlam(k,1));
    continue;
  end
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

%% Seam mask (ptt.imgCombSeam), driven by the product's own img_comb
% The band is asymmetric - the crossfade runs forward from the boundary -
% so with seam_mask_win = 1 it is [t_comb - win, t_comb + 2*win]. Coverage
% is the per-bin coherent fraction after masking, so it is exactly zero
% across the band and untouched on either side of it.
seam_lo = seam_t_comb - seam_win;
seam_hi = seam_t_comb + 2*seam_win;
band_in = @(lo,hi) Time > lo + dt/2 & Time < hi - dt/2;
band_out = @(lo,hi) (Time > lo - 6*dt & Time < lo - dt/2) | ...
  (Time > hi + dt/2 & Time < hi + 6*dt);
assert(all(all(out.coverage_blk(band_in(seam_lo,seam_hi),:) == 0)), ...
  'seam mask did not clear the crossfade band');
assert(all(all(out.coverage_blk(band_out(seam_lo,seam_hi),:) > 0.5)), ...
  'seam mask removed bins outside the declared band');

% The gap is bridged by interpolation so the inversion keeps a continuous
% chain, but the nodes that came from it must say so
assert(all(ismember(out.dlam_interpolated(:), [0;1])), ...
  'dlam_interpolated must be 0/1 for inverted blocks');
n_interp = sum(out.dlam_interpolated(:) == 1);
fprintf('Seam mask: %d bins cleared, %d of %d nodes interpolated across the gap\n', ...
  nnz(band_in(seam_lo,seam_hi)), n_interp, numel(out.dlam_interpolated));
assert(n_interp > 0, ...
  'no node was flagged as interpolated across the seam gap');

%% Regularized joint inversion (default for noisy field data)
param.fabric.inversion = 'joint';
param.fabric.out_path = 'fabric_joint';
success = fabric_task(param);
assert(success, 'fabric_task (joint mode) did not succeed');

outj = load(fullfile(outRoot, 'CSARP_fabric_joint', day_seg, ...
  sprintf('Data_%s_009.mat', day_seg)));
max_err_j = 0;
for k = 1:size(outj.dlam,1)
  if seam_gap(outj.dlam_interpolated, k), continue; end
  dmid = (outj.dlam_top_depth(k,1) + outj.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  max_err_j = max(max_err_j, abs(outj.dlam(k,1) - (lam_mid.lam(1) - lam_mid.lam(2))));
end
fprintf('Joint mode: max |error| block 1: %.3f\n', max_err_j);
assert(max_err_j < 0.05, 'joint-mode dlam deviates from truth by %.3f', max_err_j);

% Joint mode saves the per-block regularization diagnostics; stripping mode
% leaves them NaN (see ptt.invertBlocks)
assert(all(isfinite(outj.dtau_rms)) && all(isfinite(outj.reg_alpha)), ...
  'joint mode did not save dtau_rms/reg_alpha');
assert(all(ismember(outj.dlam_clipped(:), [0;1])), ...
  'joint-mode dlam_clipped must be 0/1 for inverted blocks');
assert(all(isnan(out.dtau_rms)) && all(isnan(out.reg_alpha)) ...
  && all(isnan(out.dlam_clipped(:))), ...
  'stripping mode must leave the joint-only diagnostics NaN');
fprintf('Joint diagnostics: rms %.3f ns, alpha %.3f, clipped intervals %d\n', ...
  outj.dtau_rms(1), outj.reg_alpha(1), sum(outj.dlam_clipped(:) == 1));

% Unknown solver names must fail with the documented error, not a stray
% formatting error from the message itself (mat2str rejects char in Octave)
param.fabric.inversion = 'jointt';
err = [];
try
  fabric_task(param);
catch err
end
assert(~isempty(err), 'invalid param.fabric.inversion was accepted');
assert(strcmp(err.identifier, 'ptt:invertBlocks:inversion'), ...
  'invalid inversion mode raised %s instead of ptt:invertBlocks:inversion', err.identifier);
fprintf('Invalid mode rejected: %s\n', err.message);

param.fabric.inversion = 'stripping'; % restore for subsequent sections

%% Coregistration-only mode (detected-power products: no usable phase)
param.fabric.dtau_source = 'coreg';
param.fabric.out_path = 'fabric_coreg';
success = fabric_task(param);
assert(success, 'fabric_task (coreg mode) did not succeed');

outc = load(fullfile(outRoot, 'CSARP_fabric_coreg', day_seg, ...
  sprintf('Data_%s_009.mat', day_seg)));
max_err_c = 0;
for k = 1:size(outc.dlam,1)
  if seam_gap(outc.dlam_interpolated, k), continue; end
  dmid = (outc.dlam_top_depth(k,1) + outc.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  max_err_c = max(max_err_c, abs(outc.dlam(k,1) - (lam_mid.lam(1) - lam_mid.lam(2))));
end
fprintf('Coreg-only mode: max |error| block 1: %.3f\n', max_err_c);
assert(max_err_c < 0.1, 'coreg-only dlam deviates from truth by %.3f', max_err_c);

%% Delta-k mode (split-spectrum ladder over synthetic ref/sec SLCs)
% Band-limited speckle SLCs with the true delay applied EXACTLY to the
% secondary by piecewise-constant spectral shifts (envelope + carrier,
% e^{-i 2 pi (fc + f_bb) tau}), mixed with independent noise to match the
% coherence profile. row_offset keeps its physical meaning (sec later by
% dtau), which the orientation regression uses.
fs = 1/dt;
f_bb = mod((0:Nt-1).'*fs/Nt + fs/2, fs) - fs/2;   % unshifted baseband axis
Wf = single(exp(-0.5*((abs(f_bb) - 0e6)/110e6).^8));  % band window ~+/-140 MHz
Cspec = single(complex(randn(Nt,Nx), randn(Nt,Nx))/sqrt(2)) .* Wf;
C0 = ifft(Cspec, [], 1);

tau_col = dtau_map(:,1) + timing_bias;
% Enough levels that the seam step widening the delay range does not
% coarsen the piecewise-constant quantization of the true profile
L = 32;
tau_edges = linspace(min(tau_col), max(tau_col)+eps, L+1);
sec_common = complex(zeros(Nt, Nx, 'single'));
for l = 1:L
  tau_l = (tau_edges(l) + tau_edges(l+1))/2;
  rows = tau_col >= tau_edges(l) & tau_col < tau_edges(l+1);
  if ~any(rows), continue; end
  Sl = ifft(Cspec .* single(exp(-1i*2*pi*(fc + f_bb)*tau_l)), [], 1);
  sec_common(rows,:) = Sl(rows,:);
end
clear Sl;

gam = min(max(repmat(coh_prof, 1, Nx), 0.05), 0.95);
a = single(sqrt((1 - gam)./gam));
N1 = ifft(single(complex(randn(Nt,Nx), randn(Nt,Nx))/sqrt(2)) .* Wf, [], 1);
N2 = ifft(single(complex(randn(Nt,Nx), randn(Nt,Nx))/sqrt(2)) .* Wf, [], 1);
ref = C0 + a.*N1;
sec = sec_common + a.*N2;
clear C0 sec_common N1 N2 Cspec;

% Multilooked interferogram + coherence from the SLCs (schema fields)
box = ones(9,15)/(9*15);
num = conv2(sec .* conj(ref), box, 'same');
den = sqrt(conv2(abs(ref).^2, box, 'same') .* conv2(abs(sec).^2, box, 'same'));
interferogram_mlook = single(num);
interferogram_coherence = single(abs(num)./max(den, eps));
clear num den;

in_dir2 = fullfile(outRoot, 'CSARP_polarimetric_slc', day_seg);
mkdir(in_dir2);
save('-v7', fullfile(in_dir2, sprintf('Data_%s_009.mat', day_seg)), ...
  'interferogram_mlook','interferogram_coherence','row_offset','Time', ...
  'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom', ...
  'param_records','param_polarimetric','file_type','file_version', ...
  'ref','sec');
clear ref sec interferogram_mlook interferogram_coherence;

param.fabric.dtau_source = 'deltak';
param.fabric.in_path = 'polarimetric_slc';
param.fabric.out_path = 'fabric_deltak';
success = fabric_task(param);
assert(success, 'fabric_task (deltak mode) did not succeed');

outd = load(fullfile(outRoot, 'CSARP_fabric_deltak', day_seg, ...
  sprintf('Data_%s_009.mat', day_seg)));
fprintf('Delta-k phase sign detected: %+d (expected -1)\n', outd.phase_sign);
assert(outd.phase_sign == -1, 'delta-k orientation regression failed');
max_err_d = 0;
fprintf('Delta-k block 1: interval (depth m)   true dlam   inferred\n');
for k = 1:size(outd.dlam,1)
  dmid = (outd.dlam_top_depth(k,1) + outd.dlam_bot_depth(k,1))/2;
  lam_mid = ptt.columnProfiles(parT, 1 - dmid/parT.H);
  dlam_true_k = lam_mid.lam(1) - lam_mid.lam(2);
  if seam_gap(outd.dlam_interpolated, k)
    fprintf('  %5.0f - %5.0f          %8.3f  %8.3f  (seam gap)\n', ...
      outd.dlam_top_depth(k,1), outd.dlam_bot_depth(k,1), dlam_true_k, outd.dlam(k,1));
    continue;
  end
  max_err_d = max(max_err_d, abs(outd.dlam(k,1) - dlam_true_k));
  fprintf('  %5.0f - %5.0f          %8.3f  %8.3f\n', ...
    outd.dlam_top_depth(k,1), outd.dlam_bot_depth(k,1), dlam_true_k, outd.dlam(k,1));
end
fprintf('Delta-k mode: max |error| block 1: %.3f\n', max_err_d);
assert(max_err_d < 0.06, 'delta-k dlam deviates from truth by %.3f', max_err_d);
assert(all(outd.blend_fringes == 0), ...
  'delta-k mode must not apply fringe blending');

% Delta-k resolves its ladder integers from a tau_A smoothed over the
% analysis cells and interpolated back to the full grid, so the seam band
% is widened by that reach - ((smooth-1)/2 + 1) cells, 300 ns at the
% defaults - and only for this estimator (see ptt.imgCombSeam).
dk_reach = ((5-1)/2 + 1)*100e-9;
assert(all(all(outd.coverage_blk(band_in(seam_lo-dk_reach, seam_hi+dk_reach),:) == 0)), ...
  'delta-k seam mask was not widened by the analysis-cell reach');
assert(all(all(outd.coverage_blk(band_out(seam_lo-dk_reach, seam_hi+dk_reach),:) > 0.5)), ...
  'delta-k seam mask was widened beyond the analysis-cell reach');
assert(sum(outd.dlam_interpolated(:) == 1) > 0, ...
  'no node was flagged as interpolated across the widened delta-k seam gap');
fprintf('Delta-k seam mask: %d bins cleared (%d in phase mode), %d nodes interpolated\n', ...
  nnz(band_in(seam_lo-dk_reach, seam_hi+dk_reach)), ...
  nnz(band_in(seam_lo,seam_hi)), sum(outd.dlam_interpolated(:) == 1));

fprintf('\nPASS (%.1f s)\n', toc(t0));
