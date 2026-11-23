function out = nodeFabric(S, z, dt, opts)
%NODEFABRIC One slow axis from the eigen-pair delay against azimuth, dlam(z) from the node ladder and the carrier phase, by linear ML.
%
% out = ptt.estimators.run('nodeFabric', S, z, dt, opts)
%
% The two steps of the NEGIS poster scheme ("the poster" below: orientation
% first, then one linear maximum-likelihood solve for dlam(z) from
% travel-time differences; the `nymand` method holds its original
% estimators), with the measurements this radar can carry:
%
%   ORIENTATION. The strength read from the eigen-pair carrier phase
%   (below) is not a cosine of the synthesis azimuth but a PLATEAU: a
%   pair within 45 deg of the slow axis reads the full strength (the
%   phase winds once per fringe however impure the pair), a pair within
%   45 deg of the fast axis reads its negative, and the read crosses zero
%   sharply in between (test_xcorr_traveltime verdict 5). Scanning the
%   azimuth over 0-180 deg and fitting that square wave gives the slow
%   axis directly, with no separate slow/fast step, from the same
%   measurement the strength uses. A depth-constant frame, as the poster
%   assumes. The depth swing of the 45-deg cross-pol ratio (the node
%   contrast) is kept as a diagnostic only: cross-talk puts a depth
%   modulation on the eigenframe pair too, which moves its peak.
%
%   STRENGTH. Two kinds of travel-time data enter ONE linear Gaussian
%   solve for dlam(z) (the poster's step 2):
%     - the nodes: absolute delays at the node depths, which carry the
%       fringe integer with no tracking. A node is read from POWER at
%       zero lag, so it counts the delay in half waves of the band's own
%       centre frequency (xcorrTraveltime's f_rf, measured per window),
%       NOT of the LO: n / (2 f_rf). On the accumulation radar the band
%       sits ~13 MHz below the 750 MHz LO, a 2 % difference that
%       otherwise sets the nodes against the phase;
%     - the carrier phase of the eigen-pair cross-correlation between
%       neighbouring windows, wrapped: the local slope, used only where
%       both windows are coherent and the step is small, so it never
%       crosses the odd-node coherence dips where tracking slips. Read at
%       the envelope peak, it counts the delay in fringes of the LO.
%   Both are referenced to the firn (window ref_z), which removes the
%   channel phase a co-pol gain mismatch leaves on every row.
%
%   WHAT THE PHASE CAN AND CANNOT RESOLVE. A window's phase is the
%   birefringent delay PLUS a phase the reflectors themselves put between
%   the two polarisations. That part does not accumulate: it is bounded
%   (about 3 deg of carrier on Ridge A, the same on repeat passes over
%   the same ice), so between neighbouring windows it is as large as the
%   delay step it is added to, and over a long span it is nothing. A
%   solve that treats each increment as an independent delay measurement
%   hands that phase back as 10 m stripes in dlam. The depth model
%   (opts.depth_model) therefore has two parts, both stated:
%     - dlam(z) a smooth Gaussian process about a free level (Matern 3/2,
%       sd 0.04, 500 m);
%     - a window phase psi with correlation length one window step, its
%       level set by each block's own evidence, entering the increments
%       as psi(z2) - psi(z1).
%   out.resolution_m and out.n_eff say what that leaves resolved (about
%   130 m and a dozen numbers per 1.4 km column on Ridge A), and
%   out.intervals gives the PRIOR-FREE answer at a declared span: the
%   summed increments over each interval, no model at all. The two agree
%   within their errors; where they do not, believe the intervals.
%   A change sharper than resolution_m (Ridge A's horizon near 470 m is
%   one) comes out spread over that width, and the posterior does not
%   know it: read dlam with its resolution. A ramp term for it was tried
%   and left out - in one block it is located to no better than tens of
%   metres and does not improve the levels either side; the horizon is a
%   survey-level product (scripts/figures/node_representation.py).
%   The posterior is exact for this linear model; the node residuals say
%   whether the count and the slopes agree, and the answer abstains (NaN)
%   below the first node they disagree at.
%
% Inputs
%   S     struct hh, vv, hv, vh [Nt x Nx] complex, the antenna-frame
%         channels of ONE along-track block, coregistered by a constant
%         delay per channel (never by a depth-varying field: that removes
%         the delay this measures), in the stored (engineering) sign
%   z     [Nt x 1] depth [m]
%   dt    fast-time sample interval [ns]
%   opts  (all optional)
%     .calibrate    (true)    put the four channels on one complex gain
%                   scale first (ptt.channelGains: reciprocity + firn).
%                   Uncalibrated, the synthesised eigen pair is impure and
%                   the phase increments carry a strength bias that
%                   alternates with the node ladder.
%     .gains        struct gr, gt to apply instead of estimating them
%     .firn_m       ([30 90]) the firn window for the co-pol gain
%     .psi0         slow-axis azimuth in the ANTENNA frame [rad]; when
%                   given the scan is skipped
%     .psi_step_deg (5)       orientation scan step
%     .z_band       ([300 1450]) depths the orientation contrast is read over
%     .fc           (750e6)   LO for the carrier phase (xcorrTraveltime)
%     .node_fc      ([])      frequency the nodes count in [Hz]; empty =
%                             measured from the data (windows.f_rf)
%     .conj         (true)    stored sign (xcorrTraveltime)
%     .win_m        (50)      xcorr window; step_m (25)
%     .ref_z        ([30 90]) the firn reference window [m]
%     .gamma_min    (0.2)     increments need both windows at or above this
%     .inc_max_frac (0.25)    increments above this fraction of a fringe are dropped
%     .sigma_scale  (1.3)     inflation of the xcorr sigma (measured optimism)
%     .sigma_floor  (0.005)   [ns] floor on increment sigma
%     .node         struct    passed to ptt.crosspolNodes
%     .depth_model  struct, fields override the defaults:
%                     kernel 'matern32' | 'exp'; sd (0.04); len_m (500);
%                     mean (0); level_sd (0.2, a free level; 0 = none);
%                     psi_sd_ns ('auto' = the block's maximum evidence, or
%                     a number, 0 = no window phase); psi_len_m (step_m)
%                   or the name 'poster' for the poster's own model
%     .prior        [mean sd len_m] the poster's model: an exponential
%                   prior about `mean`, independent increments, no window
%                   phase. With node_fc = fc this reproduces the products
%                   made before the depth model existed.
%     .model_dz     (10)      [m] dlam model grid
%     .interval_m   (125)     [m] span of the prior-free interval means
%     .node_tol_hf  (0.4)     a candidate extremum is a node only if the
%                             phase-integrated delay puts a node of its
%                             parity within this many half waves of it
%     .node_tol_frac (0.25)   node residual above this fraction of a fringe
%                             is a disagreement: abstain below it
%
% Output out
%   .gains, .gain_info  what was applied (ptt.channelGains)
%   .psi0, .theta_scan: psis_deg, strength_read (the band's strength as
%           read at each azimuth), template_score, theta_slow_deg,
%           plateau_ratio (plateau separation over its scatter: the
%           orientation's quality), contrast_45_db and contrast_0_db
%           (the node-contrast diagnostic 45 deg off the axis and on it)
%   .windows  the xcorrTraveltime output at psi0 (eigen pair)
%   .ratio45, .ratio0 [Nt x 1] cross/co ratio in dB 45 deg off and in frame
%   .nodes    ptt.crosspolNodes output, plus n, fc_hz, tau_ns, resid_ns, parity_ok
%   .fine_ref [ns] the firn phase reference
%   .zm, .dlam, .dlam_raw, .sigma_dlam, .chi2_dof, .n_inc, .n_abs
%   .resolution    diagonal of the resolution matrix
%   .resolution_m  width [m] holding the central 68 % of each cell's
%                  averaging kernel: what one plotted value averages over
%   .n_eff         trace of the resolution matrix: how many independent
%                  numbers the data put into the profile
%   .depth_model   the model used, psi_sd_ns as estimated
%   .psi           [Nw x 1] posterior window phase [ns] (NaN where a window
%                  enters no increment)
%   .intervals     z1, z2, dlam, sigma, coverage: the summed increments
%                  over each interval_m span, no prior; NaN below 60 %
%                  coverage
%   .fit           resid_acf1: lag-one autocorrelation of the increment
%                  residuals. About -0.3 when what is left is a bounded
%                  window phase; near zero or positive means structure
%                  that ACCUMULATES is being smoothed over.
%   .log_evidence  of the final solve
%   .abstain_below_m  depth below which dlam is NaN (Inf when none)
%   .k_ns_per_m
%   .n_traces_dropped traces left out because a channel had non-finite samples
if nargin < 4, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
conj_ = H_opt(opts, 'conj', true);
step_deg = H_opt(opts, 'psi_step_deg', 5);
z_band = H_opt(opts, 'z_band', [300 1450]);
win_m = H_opt(opts, 'win_m', 50);
step_m = H_opt(opts, 'step_m', 25);
ref_z = H_opt(opts, 'ref_z', [30 90]);
gamma_min = H_opt(opts, 'gamma_min', 0.2);
inc_max_frac = H_opt(opts, 'inc_max_frac', 0.25);
sigma_scale = H_opt(opts, 'sigma_scale', 1.3);
sigma_floor = H_opt(opts, 'sigma_floor', 0.005);
node_opts = H_opt(opts, 'node', struct());
model = H_depth_model(opts, step_m);
model_dz = H_opt(opts, 'model_dz', 10);
interval_m = H_opt(opts, 'interval_m', 125);
node_tol = H_opt(opts, 'node_tol_frac', 0.25);
node_tol_hf = H_opt(opts, 'node_tol_hf', 0.4);
node_fc_fixed = H_opt(opts, 'node_fc', []);
z = double(z(:));
dz = median(diff(z));
fringe = 1e9 / fc;
% ---- a trace any channel lacks is no use to a four-channel synthesis: drop
% it here, once, or one gap in the cache turns every trace-wise mean below
% into NaN (two Ridge A frames came back with no gains and no nodes from a
% single all-NaN trace). A block missing most of its traces is no block.
keep_x = all(isfinite(S.hh), 1) & all(isfinite(S.vv), 1) & all(isfinite(S.hv), 1) & all(isfinite(S.vh), 1);
n_dropped = nnz(~keep_x);
if n_dropped > 0
  if nnz(keep_x) < max(8, 0.5 * numel(keep_x))
    error('ptt:nodeFabric:traces', '%d of %d traces have non-finite samples', n_dropped, numel(keep_x));
  end
  for ch = {'hh', 'vv', 'hv', 'vh'}, S.(ch{1}) = S.(ch{1})(:, keep_x); end
end
C = ptt.constants();
k = C.deps / (sqrt(C.eps_bar) * C.c);                 % ns per m per unit dlam, two-way
xo = struct('fc', fc, 'conj', conj_, 'win', max(8, round(win_m / dz)), 'step', max(1, round(step_m / dz)));

% ---- one gain scale for the four channels
gains = struct('gr', 1, 'gt', 1); gain_info = struct();
if isfield(opts, 'gains') && ~isempty(opts.gains)
  gains = opts.gains;
elseif H_opt(opts, 'calibrate', true)
  [gains, gain_info] = ptt.channelGains(S, z, struct('firn_m', H_opt(opts, 'firn_m', [30 90])));
end
if ~(gains.gr == 1 && gains.gt == 1)
  S.vv = S.vv / (gains.gr * gains.gt);
  S.hv = S.hv / gains.gr;
  S.vh = S.vh / gains.gt;
end

% ---- orientation
scan = struct('psis_deg', [], 'strength_read', [], 'template_score', [], ...
  'theta_slow_deg', NaN, 'plateau_ratio', NaN, 'contrast_45_db', NaN, 'contrast_0_db', NaN);
if isfield(opts, 'psi0') && ~isempty(opts.psi0)
  psi0 = opts.psi0;
else
  psis = 0:step_deg:180 - step_deg;
  sread = nan(size(psis));
  band = z >= z_band(1) & z <= z_band(2);
  rows = find(z >= z_band(1) - 100 & z <= z_band(2) + 100);     % the band plus a window either side
  Sb = struct('hh', S.hh(rows, :), 'vv', S.vv(rows, :), 'hv', S.hv(rows, :), 'vh', S.vh(rows, :));
  for i = 1:numel(psis)
    [a, b] = H_copair(Sb, deg2rad(psis(i)));
    w = ptt.xcorrTraveltime(a, b, dt, xo);
    [inc, ~, zi, z1, z2] = H_increments(w, z(rows), fringe, gamma_min, inc_max_frac, 1, 0);
    in = zi >= z_band(1) & zi <= z_band(2);
    if any(in), sread(i) = sum(inc(in)) / (k * sum(z2(in) - z1(in))); end
  end
  % the read is a plateau: +S within 45 deg of the slow axis, -S within
  % 45 deg of the fast one, crossing zero in between. Fit that square wave
  % for the axis directly - the fit's sign IS the slow/fast decision
  th = 0:0.5:179.5;
  score = nan(size(th));
  ok = isfinite(sread);
  for j = 1:numel(th)
    d = abs(mod(psis - th(j) + 90, 180) - 90);
    score(j) = sum(sread(ok) .* (2 * (d(ok) < 45) - 1));
  end
  [~, jb] = max(score);
  theta_slow = th(jb);
  % refine from the two zero crossings, interpolated between the samples
  % that straddle them: the fast edge at theta - 45 (read turning
  % positive) and theta + 45 (turning negative)
  est = [];
  for side = [-1, 1]
    zc0 = theta_slow + side * 45;
    for i = 1:numel(psis)
      i2 = mod(i, numel(psis)) + 1;
      p1 = psis(i); p2 = psis(i) + step_deg;
      if ~(isfinite(sread(i)) && isfinite(sread(i2))), continue; end
      if abs(mod(p1 - zc0 + 90, 180) - 90) > step_deg, continue; end
      if sread(i) * sread(i2) < 0 && sign(sread(i2) - sread(i)) == -side
        x0 = p1 + (p2 - p1) * sread(i) / (sread(i) - sread(i2));
        est(end+1) = x0 - side * 45; %#ok<AGROW>
      end
    end
  end
  if ~isempty(est)
    theta_slow = mod(rad2deg(angle(mean(exp(2i * deg2rad(est))))) / 2, 180);
  end
  d = abs(mod(psis - theta_slow + 90, 180) - 90);
  s_in = sread(ok & d < 40); s_out = sread(ok & d > 50);
  plateau = (mean(s_in) - mean(s_out)) / max(2 * std([s_in, -s_out]), eps);
  psi0 = deg2rad(theta_slow);
  % the node-contrast diagnostic: the depth swing of the 45-deg cross-pol
  % ratio at the slow axis's 45-deg pair against the eigenframe pair's
  % own (ideally none; cross-talk puts some there)
  con = [H_contrast(H_ratio(S, psi0 + pi / 4), z, band, dz), H_contrast(H_ratio(S, psi0), z, band, dz)];
  scan = struct('psis_deg', psis, 'strength_read', sread, 'template_score', score, ...
    'theta_slow_deg', theta_slow, 'plateau_ratio', plateau, 'contrast_45_db', con(1), 'contrast_0_db', con(2));
end

% ---- measurements at psi0
T0 = ptt.rotatePolarization(S, psi0);
w = ptt.xcorrTraveltime(T0.hh, T0.vv, dt, xo);
zw = z(w.idx);
ratio45 = H_ratio(S, psi0 + pi / 4);
ratio0 = H_ratio(S, psi0);
N = ptt.crosspolNodes(ratio45, z, node_opts);
% the firn phase reference: a co-pol gain phase sits on every row, and the
% column above ref_z has no birefringence to speak of
inref = zw >= ref_z(1) & zw <= ref_z(2) & isfinite(w.fine);
if ~any(inref), inref = (1:numel(zw)).' <= 2; end
fine_ref = angle(mean(exp(2i * pi * w.fine(inref) / fringe))) / (2 * pi) * fringe;
zref = mean(zw(inref));
% the frequency a node counts the delay in: the band's own centre where the
% node is, from the windows around it (a running median over ~200 m; one
% window's centroid wanders a few MHz). Not the LO - see the header.
if ~isempty(node_fc_fixed)
  f_node = @(zz) node_fc_fixed * ones(size(zz));
else
  fw = w.f_rf;
  okf = isfinite(fw) & w.gamma >= gamma_min;
  if nnz(okf) >= 3
    fs = movmedian(fw(okf), 9);
    zf = zw(okf);
    f_node = @(zz) interp1(zf, fs, min(max(zz, zf(1)), zf(end)));
  else
    f_node = @(zz) fc * ones(size(zz));
  end
end

% ---- data rows
[inc, sig_inc, zi, z1, z2, iw] = H_increments(w, z, fringe, gamma_min, inc_max_frac, sigma_scale, sigma_floor);
% ---- model
zm = (0:model_dz:max(z)).';
zm = zm(zm < max(z));
Nm = numel(zm);
A_abs = @(zz) k * max(0, min(model_dz, reshape(zz, [], 1) - zm.'));   % cumulative thickness above each depth
m_prior = model.mean * ones(Nm, 1);
Cm = H_kernel(zm, model);
G_inc = A_abs(z2) - A_abs(z1);
G0 = [A_abs(zref); G_inc]; d0 = [0; inc];
% the window phase's level, from this block's own evidence: first on the
% phase rows alone, which is all there is until the nodes are counted
psi_auto = ischar(model.psi_sd_ns) || isstring(model.psi_sd_ns);
psi_curve = [];
if psi_auto
  model.psi_sd_ns = H_psi_auto(G0, d0, 0.02, sig_inc, iw, zw, model.psi_len_m, Cm, m_prior);
end
% ---- which candidates are nodes, and their count, from the phase. The
% increments integrated from the surface give a delay at every candidate;
% a maximum is a node only where that delay is an odd number of half
% waves, a minimum an even one, the count strictly increasing down the
% ladder. Each candidate is judged against the previous ACCEPTED node, so
% a coherence dip the increments had to bridge cannot throw off every
% node below it. Speckle wiggles on a broad extremum fail the parity or
% the spacing and are set aside (out.nodes.rejected).
n_source = 'none';
cand = N;
nn = nan(size(N.z)); acc = false(size(N.z)); hf = nan(size(N.z));
fn_c = f_node(N.z);
if ~isempty(N.z)
  m0 = H_solve(G0, d0, H_cd(0.02, sig_inc, iw, zw, model.psi_sd_ns, model.psi_len_m), Cm, m_prior);
  hf = 2 * (A_abs(N.z) * m0) .* fn_c * 1e-9;            % half waves at each candidate
  prev_n = 0; prev_hf = 0;
  for j = 1:numel(N.z)
    odd = N.type(j) > 0;
    c = prev_n + (1:8);
    c = c(mod(c, 2) == odd);
    [off, i] = min(abs(c - prev_n - (hf(j) - prev_hf)));
    if off <= node_tol_hf
      nn(j) = c(i); acc(j) = true;
      prev_n = nn(j); prev_hf = hf(j);
    end
  end
  if any(acc), n_source = 'phase'; end
end
N = struct('z', cand.z(acc), 'zsig', cand.zsig(acc), 'type', cand.type(acc), 'n', nn(acc), ...
  'prom_db', cand.prom_db(acc), 'hf_phase', hf(acc), 'fc_hz', fn_c(acc), 'first_ok', any(acc), ...
  'n_source', n_source, 'smoothed', cand.smoothed, ...
  'rejected', struct('z', cand.z(~acc), 'type', cand.type(~acc), 'prom_db', cand.prom_db(~acc)));
good = isfinite(N.n);
% a node's delay: n half waves of the band centre, less the firn phase
% (the same phase, so the same number of cycles, at the node's frequency)
tau_node = N.n ./ (2 * N.fc_hz * 1e-9) - fine_ref * fc ./ N.fc_hz;
sig_node = max(0.05, k * 0.06 * N.zsig);                 % the depth uncertainty at a nominal strength
% ---- solve with nodes and increments
G = [A_abs(zref); A_abs(N.z(good)); G_inc];
d = [0; tau_node(good); inc];
sd_abs = [0.02; sig_node(good)];
% the window phase level again, now that every row is in
if psi_auto
  [model.psi_sd_ns, psi_curve] = H_psi_auto(G, d, sd_abs, sig_inc, iw, zw, model.psi_len_m, Cm, m_prior);
end
Cd = H_cd(sd_abs, sig_inc, iw, zw, model.psi_sd_ns, model.psi_len_m);
[m, Cpost, Rm, res, chi2, dof, logev, alpha] = H_solve(G, d, Cd, Cm, m_prior);
n_abs = 1 + nnz(good); n_inc = numel(inc);
% the window phase the solve set aside, and what is left on the increments
psi = nan(numel(zw), 1);
acf1 = NaN;
if n_inc > 0
  r_inc = res(n_abs + 1:end);
  if model.psi_sd_ns > 0
    [~, Dw, Rp] = H_cd(sd_abs, sig_inc, iw, zw, model.psi_sd_ns, model.psi_len_m);
    p_hat = model.psi_sd_ns ^ 2 * (Rp * (Dw.' * alpha(n_abs + 1:end)));
    used = false(numel(zw), 1); used(iw) = true; used(iw + 1) = true;
    psi(used) = p_hat(used);
  end
  nb = find(diff(iw) == 1);                              % neighbouring increments
  if numel(nb) >= 5 && any(r_inc)
    acf1 = mean(r_inc(nb) .* r_inc(nb + 1)) / mean(r_inc .^ 2);
  end
end
% node residuals and abstention
nn = N.n;
node_res = nan(size(nn));
node_res(good) = res(2:1 + nnz(good));
bad = good & abs(node_res) > node_tol * fringe;
abstain = Inf;
if any(bad), abstain = min(N.z(bad)); end
if n_inc == 0 && ~any(good), abstain = 0; end     % no datum at all: the prior is not an answer
dlam = m; dlam(zm >= abstain) = NaN;
% parity of the carrier phase at the nodes: a node is where the pair's own
% phase, firn offset and all, is n half turns, so the carrier phase delay
% there should be n half waves to within a sixth of a fringe (unreliable at
% odd nodes, where the pair decorrelates; reported, not used)
par = nan(size(nn));
for i = find(good).'
  f = interp1(zw(isfinite(w.fine)), unwrap_free(w.fine(isfinite(w.fine)), fringe), N.z(i), 'nearest', NaN);
  if isfinite(f)
    ph = mod(f - nn(i) / (2 * N.fc_hz(i) * 1e-9) + fringe / 2, fringe) - fringe / 2;
    par(i) = abs(ph) < fringe / 6;
  end
end
N.tau_ns = tau_node; N.sigma_ns = sig_node; N.resid_ns = node_res; N.parity_ok = par;
% what one plotted value averages over
res_m = nan(Nm, 1);
for i = 1:Nm
  a = abs(Rm(i, :));
  if sum(a) > 0
    ca = cumsum(a) / sum(a);
    [cu, iu] = unique(ca);
    if numel(cu) >= 2, res_m(i) = diff(interp1(cu, zm(iu), [0.16 0.84], 'linear', 'extrap')); end
  end
end

out = struct('gains', gains, 'gain_info', gain_info, 'psi0', psi0, 'theta_scan', scan, 'windows', w, 'ratio45', ratio45, 'ratio0', ratio0, ...
  'nodes', N, 'fine_ref', fine_ref, 'zm', zm, 'dlam', dlam, 'dlam_raw', m, ...
  'sigma_dlam', sqrt(max(0, diag(Cpost))), 'resolution', diag(Rm), 'resolution_m', res_m, 'n_eff', trace(Rm), ...
  'chi2_dof', chi2 / dof, 'n_inc', n_inc, 'n_abs', n_abs, 'abstain_below_m', abstain, 'k_ns_per_m', k, ...
  'inc', struct('z1', z1, 'z2', z2, 'dtau', inc, 'sigma', sig_inc, 'zc', zi, 'iw', iw), 'fringe_ns', fringe, ...
  'depth_model', model, 'psi', psi, 'psi_evidence', psi_curve, ...
  'intervals', H_intervals(inc, z1, z2, zi, zw, iw, sig_inc, model, k, interval_m), ...
  'fit', struct('resid_acf1', acf1), 'log_evidence', logev, ...
  'n_traces_dropped', n_dropped);
end

% ---------------------------------------------------------------------
function M = H_depth_model(opts, step_m)
% the depth model in force: the defaults, the poster's model if the legacy
% opts.prior (or the name 'poster') asks for it, then any fields given
M = struct('kernel', 'matern32', 'sd', 0.04, 'len_m', 500, 'mean', 0, 'level_sd', 0.2, ...
  'psi_sd_ns', 'auto', 'psi_len_m', step_m);
dm = H_opt(opts, 'depth_model', []);
pr = H_opt(opts, 'prior', []);
if (ischar(dm) || isstring(dm)) && strcmpi(dm, 'poster'), pr = [0 0.05 150]; dm = []; end
if (ischar(dm) || isstring(dm)) && ~strcmpi(dm, 'smooth')
  error('ptt:nodeFabric:model', 'depth_model must be a struct, ''smooth'' or ''poster''');
end
if ~isempty(pr)
  M = struct('kernel', 'exp', 'sd', pr(2), 'len_m', pr(3), 'mean', pr(1), 'level_sd', 0, ...
    'psi_sd_ns', 0, 'psi_len_m', step_m);
end
if isstruct(dm)
  fn = fieldnames(dm);
  for i = 1:numel(fn)
    if ~isfield(M, fn{i}), error('ptt:nodeFabric:model', 'unknown depth_model field %s', fn{i}); end
    M.(fn{i}) = dm.(fn{i});
  end
end
end

function Cm = H_kernel(zm, M)
% prior covariance of dlam on the model grid
D = abs(zm - zm.');
switch lower(M.kernel)
  case 'exp'
    Cm = M.sd ^ 2 * exp(-D / M.len_m);
  case 'matern32'
    a = sqrt(3) * D / M.len_m;
    Cm = M.sd ^ 2 * (1 + a) .* exp(-a);
  otherwise
    error('ptt:nodeFabric:model', 'unknown kernel %s', M.kernel);
end
Cm = Cm + 1e-12 * M.sd ^ 2 * eye(numel(zm)) + M.level_sd ^ 2;
end

function [Cd, D, Rp] = H_cd(sd_abs, sd_inc, iw, zw, psi_sd, psi_len)
% data covariance: the absolute rows (firn anchor, nodes) independent; the
% increments their own noise plus the window phase differenced,
% psi(z2) - psi(z1), which ties neighbouring increments together
n_inc = numel(sd_inc);
Ci = diag(sd_inc(:) .^ 2);
D = []; Rp = [];
if psi_sd > 0 && n_inc > 0
  nw = numel(zw);
  D = sparse([1:n_inc, 1:n_inc], [iw(:).', iw(:).' + 1], [-ones(1, n_inc), ones(1, n_inc)], n_inc, nw);
  Rp = exp(-abs(zw(:) - zw(:).') / max(psi_len, eps));
  Ci = Ci + psi_sd ^ 2 * full(D * Rp * D.');
end
Cd = blkdiag(diag(sd_abs(:) .^ 2), Ci);
end

function [m, Cpost, Rm, res, chi2, dof, logev, alpha] = H_solve(G, d, Cd, Cm, m_prior)
% the linear Gaussian solve in data space,
%   m = m_prior + Cm G' (G Cm G' + Cd)^-1 (d - G m_prior),
% which needs no inverse of the prior (a smooth kernel has none worth the
% name) and gives the evidence for free
Sd = G * Cm * G.' + Cd;
Sd = (Sd + Sd.') / 2;
Lc = chol(Sd, 'lower');
r0 = d - G * m_prior;
alpha = Lc.' \ (Lc \ r0);
K = Cm * G.';
m = m_prior + K * alpha;
if nargout < 2, return; end
Rm = K * (Lc.' \ (Lc \ G));
Cpost = Cm - K * (Lc.' \ (Lc \ K.'));
res = d - G * m;
chi2 = res.' * (Cd \ res);
dof = max(1, numel(d) - trace(Rm));
logev = -0.5 * (r0.' * alpha) - sum(log(diag(Lc))) - 0.5 * numel(d) * log(2 * pi);
end

function le = H_logev(G, d, Cd, Cm, m_prior)
Sd = G * Cm * G.' + Cd;
Lc = chol((Sd + Sd.') / 2, 'lower');
y = Lc \ (d - G * m_prior);
le = -0.5 * (y.' * y) - sum(log(diag(Lc))) - 0.5 * numel(d) * log(2 * pi);
end

function [psi_sd, curve] = H_psi_auto(G, d, sd_abs, sd_inc, iw, zw, psi_len, Cm, m_prior)
% the window phase level that maximises this block's evidence, on a log
% grid with a parabola through the top three; zero when no level beats none
grid = logspace(log10(0.002), log10(0.04), 14);
le = nan(size(grid));
for i = 1:numel(grid)
  le(i) = H_logev(G, d, H_cd(sd_abs, sd_inc, iw, zw, grid(i), psi_len), Cm, m_prior);
end
le0 = H_logev(G, d, H_cd(sd_abs, sd_inc, iw, zw, 0, psi_len), Cm, m_prior);
curve = struct('psi_sd_ns', [0, grid], 'log_evidence', [le0, le]);
[le_max, i] = max(le);
if numel(sd_inc) < 5 || le0 >= le_max
  psi_sd = 0;
  return
end
psi_sd = grid(i);
if i > 1 && i < numel(grid)
  x = log(grid(i-1:i+1)); p = polyfit(x - x(2), le(i-1:i+1), 2);
  if p(1) < 0, psi_sd = exp(x(2) - p(2) / (2 * p(1))); end
end
end

function I = H_intervals(inc, z1, z2, zi, zw, iw, sig_inc, model, k, span)
% dlam over fixed spans from the summed increments alone: no prior, no
% nodes. The error is the increments' own covariance summed, so the window
% phase enters once at each end however long the span; an interval less
% than 60 % covered by increments is left NaN rather than filled in.
I = struct('z1', [], 'z2', [], 'dlam', [], 'sigma', [], 'coverage', []);
if isempty(inc) || ~(span > 0), return; end
e = (span * ceil(zw(1) / span):span:zw(end)).';
if numel(e) < 2, return; end
ps = model.psi_sd_ns; if ~isnumeric(ps), ps = 0; end
Ci = H_cd(zeros(0, 1), sig_inc, iw, zw, ps, model.psi_len_m);
n = numel(e) - 1;
I.z1 = e(1:n); I.z2 = e(2:n + 1);
I.dlam = nan(n, 1); I.sigma = nan(n, 1); I.coverage = zeros(n, 1);
for j = 1:n
  in = zi >= e(j) & zi < e(j + 1);
  th = sum(z2(in) - z1(in));
  I.coverage(j) = th / span;
  if I.coverage(j) >= 0.6
    I.dlam(j) = sum(inc(in)) / (k * th);
    I.sigma(j) = sqrt(sum(Ci(in, in), 'all')) / (k * th);
  end
end
end

function [a, b] = H_copair(S, psi)
% the co-polarised pair at azimuth psi only - rotatePolarization's T_hh
% and T_vv - without forming the cross-polarised channels the scan does
% not need
c = cos(psi); s_ = sin(psi);
x = double(S.hv) + double(S.vh);
a = c ^ 2 * double(S.hh) + s_ ^ 2 * double(S.vv) + c * s_ * x;
b = s_ ^ 2 * double(S.hh) + c ^ 2 * double(S.vv) - c * s_ * x;
end

function r = H_ratio(S, psi)
% cross/co power ratio [dB] of the pair synthesised at azimuth psi, pooled
% over traces
T = ptt.rotatePolarization(S, psi);
px = mean(abs(T.hv) .^ 2 + abs(T.vh) .^ 2, 2) / 2;
pc = mean(abs(T.hh) .^ 2 + abs(T.vv) .^ 2, 2) / 2;
r = 10 * log10(max(px, realmin) ./ max(pc, realmin));
end

function c = H_contrast(r, z, band, dz)
% rms of the smoothed ratio about a quadratic trend over the band
rs = movmean(r, max(1, round(20 / dz)));
x = z(band); y = rs(band);
ok = isfinite(y);
if nnz(ok) < 10, c = NaN; return; end
p = polyfit((x(ok) - mean(x(ok))) / 1000, y(ok), 2);
c = rms(y(ok) - polyval(p, (x(ok) - mean(x(ok))) / 1000));
end

function [inc, sig, zc, z1, z2, iw] = H_increments(w, z, fringe, gamma_min, inc_max_frac, sigma_scale, sigma_floor)
% wrapped carrier-phase increments between neighbouring windows; iw is the
% index of each increment's upper window (the lower one is iw + 1)
zw = z(w.idx);
f = w.fine; g = w.gamma; s = w.sigma;
inc = zeros(0, 1); sig = zeros(0, 1); zc = zeros(0, 1); z1 = zeros(0, 1); z2 = zeros(0, 1); iw = zeros(0, 1);
for i = 1:numel(f) - 1
  if ~(g(i) >= gamma_min && g(i+1) >= gamma_min && isfinite(f(i)) && isfinite(f(i+1))), continue; end
  d = mod(f(i+1) - f(i) + fringe / 2, fringe) - fringe / 2;
  if abs(d) > inc_max_frac * fringe, continue; end
  inc(end+1, 1) = d; %#ok<AGROW>
  sig(end+1, 1) = max(sigma_floor, sigma_scale * hypot(s(i), s(i+1))); %#ok<AGROW>
  zc(end+1, 1) = (zw(i) + zw(i+1)) / 2; %#ok<AGROW>
  z1(end+1, 1) = zw(i); z2(end+1, 1) = zw(i+1); %#ok<AGROW>
  iw(end+1, 1) = i; %#ok<AGROW>
end
end

function f = unwrap_free(f, fringe)
% identity: the fine delay is used wrapped; kept as a named step so the
% parity check reads as what it is
f = mod(f + fringe / 2, fringe) - fringe / 2;
end

function v = H_opt(o, name, default)
if isfield(o, name) && ~isempty(o.(name)), v = o.(name); else, v = default; end
end
