%TEST_NODE_FABRIC The node method: one slow axis per column, dlam(z) from nodes plus carrier phase.
%
% nodeFabric is the NEGIS poster's two-step scheme (orientation first, then
% a linear solve for dlam(z) from travel-time differences) with the
% measurements this
% radar can carry: the half-wave node ladder of the 45-deg cross-pol power
% (the fringe integer, no tracking) and the eigen-pair carrier phase
% between windows (the local slope), in one linear Gaussian solve whose
% depth model is a smooth dlam(z) plus a bounded window phase; orientation
% from the azimuth scan of the strength read. Tested on four-channel
% synthetic columns (synth_quadpol_column) with speckle, depolarisation,
% an off-centre band, uncalibrated channel gains and a layer-locked
% reflection phase.
%
% The verdicts, in the order that matters:
%
%   1. THE DELAY SCALE. The ns-per-metre constant is Rathmann's nadir
%      value, deps / (sqrt(eps_bar) c) = 0.06369 ns per metre of ice per
%      unit dlam, two-way.
%   2. ORIENTATION AND SLOW AXIS, CLEAN COLUMN. Axes at 35 and at 125 deg
%      from the antenna frame (the 90-deg alias pair) both come back
%      within a few degrees WITH the slow one named.
%   3. NODES AND STRENGTH, CLEAN COLUMN. Node depths within 15 m of where
%      the true delay below the firn reference passes n half waves of the
%      band's centre less the firn phase; dlam(z) within its posterior, no
%      node disagreement, no abstention.
%   4. UNCALIBRATED CHANNELS. Receive and transmit gains with a 45-deg
%      co-pol phase and the -111-deg HV/VH phase this system shows, with
%      depolarisation 0.5 and noise: orientation still within 5 deg, the
%      nodes within 30 m of the ladder (the weak contrast pulls the
%      extrema ~20 m shallow), dlam recovered. This is the verdict that
%      matters for real data.
%   5. WEAK FABRIC, NO NODES. A column whose delay never reaches half a
%      fringe has no node ladder; dlam comes from the increments alone,
%      first_ok reports the empty ladder, nothing abstains.
%   6. A MISSING TRACE. Two Ridge A frames came back with no gains and no
%      nodes from ONE all-NaN trace in their coregistered cache: every
%      trace-wise mean below it turned NaN. Such traces are dropped at
%      the door and the answer must not move.
%   7. THE NODES COUNT IN THE BAND'S FREQUENCY, NOT THE LO'S. A node is
%      read from power at zero lag, so it sits where the delay is n half
%      waves of the band's centre (13 MHz below the LO here, as on the
%      radar). On a nine-node column the true delay at the found nodes
%      matches n / (2 f_band) and misses n / (2 fc) by a few hundredths of
%      a nanosecond per node, growing down the ladder; the frequency is
%      measured from the data to within a MHz and tracks the band when it
%      is moved.
%   8. A REFLECTION PHASE IS NOT FABRIC. With a layer-locked phase between
%      the two polarisations (4.5 deg, 25 m: what Ridge A shows), the
%      poster's model - every increment an independent delay measurement
%      - hands it back as 10 m stripes as large as on the real sections.
%      The depth model returns the smooth truth inside its posterior,
%      finds the phase level from the block's own evidence, and leaves
%      residuals with the negative lag-one correlation of a bounded
%      phase. The prior-free interval means are right within the errors
%      they state.
%   9. WINDOWS. The acceptance test of 27 Sep: doubling the correlation
%      window (50 -> 100 m, step 25 -> 50 m) must not move the answer.
%      Under the depth model it moves by less than the posterior; under
%      the poster's model the 10 m profile moves by its whole stripe
%      amplitude.
%  10. THE GRID. dlam is solved on a 10 m grid. On a 20 m, a 5 m and a
%      2.5 m grid the profile is the same to a small fraction of its
%      posterior, with the same number of independent values: the grid
%      sets where the answer is written, not what it resolves.
%
% Run: matlab -batch "run('opr_fabric/test/test_node_fabric.m')"
clear; t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
addpath(thisDir);
rng(11);

fs = 300e6; dt = 1e9 / fs; fc = 750e6; f_off = -13e6;
C = ptt.constants(); n_ice = sqrt(C.eps_bar);
dz = dt * C.c / (2 * n_ice);
Nt = 5000; z = (0:Nt-1).' * dz;
dl_true = 0.012 + 0.058 * (z / 1400) .^ 1.5 + 0.015 * exp(-((z - 650) / 60) .^ 2);
fringe = 1e9 / fc;
half_wave = 1e9 / (2 * (fc + f_off));         % ns between nodes: the band's centre, not the LO
p = struct('fc', fc, 'f_off', f_off, 'Nx', 300, 'depol', 0, 'snr_db', 15, 'tau_x', 3);
opts = struct('fc', fc, 'psi_step_deg', 5);
% the poster's model as it ran before the depth model existed
poster = struct('prior', [0 0.05 150], 'node_fc', fc);

%% ---- verdict 1: the delay scale
[S, tau_true] = synth_quadpol_column(z, dl_true, deg2rad(35), dt, p);
o = ptt.nodeFabric(S, z, dt, opts);
k_ref = C.deps / (sqrt(C.eps_bar) * C.c);
fprintf('1. k: nodeFabric %.6f ns/m, deps / (sqrt(eps_bar) c) %.6f ns/m\n', o.k_ns_per_m, k_ref);
ok1 = abs(o.k_ns_per_m - k_ref) < 1e-12 && abs(k_ref - 0.063688) < 1e-6;
fprintf('   the delay scale is Rathmann''s nadir constant:             %s\n', H_tick(ok1));

%% ---- verdict 2: orientation and slow axis, both aliases
th35 = rad2deg(o.psi0);
S2 = synth_quadpol_column(z, dl_true, deg2rad(125), dt, p);
o2 = ptt.nodeFabric(S2, z, dt, opts);
th125 = rad2deg(o2.psi0);
e35 = abs(angle(exp(2i * deg2rad(th35 - 35)))) / 2 * 180 / pi;
e125 = abs(angle(exp(2i * deg2rad(th125 - 125)))) / 2 * 180 / pi;
fprintf('\n2. slow axis: truth 35 -> %.1f (err %.1f), truth 125 -> %.1f (err %.1f); plateau ratio %.1f / %.1f\n', ...
  th35, e35, th125, e125, o.theta_scan.plateau_ratio, o2.theta_scan.plateau_ratio);
ok2 = e35 < 3 && e125 < 3;
fprintf('   orientation and slow axis, both aliases:                 %s\n', H_tick(ok2));

%% ---- verdict 3: nodes and strength, clean column
zn_true = H_ladder(o, z, tau_true, half_wave);
N = o.nodes;
fprintf('\n3. nodes found %d (first_ok %d), true %d; depths:\n', numel(N.z), N.first_ok, numel(zn_true));
fprintf('   found %s\n   true  %s\n', sprintf('%5.0f ', N.z), sprintf('%5.0f ', zn_true));
[zerr, extra_ok] = H_node_check(N, zn_true, fringe);
dl_t = interp1(z, dl_true, o.zm + 5, 'linear', 'extrap');   % cell centres
ok = isfinite(o.dlam) & o.zm > 100;
zsc = (o.dlam(ok) - dl_t(ok)) ./ o.sigma_dlam(ok);
fprintf('   node depth err max %.0f m, node resid max %.3f ns, dlam median |err| %.4f, z rms %.2f, chi2/dof %.2f, abstain below %g m, inc %d\n', ...
  zerr, max(abs(N.resid_ns)), median(abs(o.dlam(ok) - dl_t(ok))), rms(zsc), o.chi2_dof, o.abstain_below_m, o.n_inc);
ok3 = N.first_ok && zerr < 15 && extra_ok && median(abs(o.dlam(ok) - dl_t(ok))) < 0.004 && ...
  rms(zsc) < 2 && isinf(o.abstain_below_m);
fprintf('   nodes and strength on the clean column:                  %s\n', H_tick(ok3));

%% ---- verdict 4: uncalibrated channels, depolarised
gr = 0.8 * exp(-1i * deg2rad(33)); gt = 0.9 * exp(1i * deg2rad(78));   % arg(gr gt) = 45, arg(gr/gt) = -111
p4 = p; p4.depol = 0.5; p4.gr = gr; p4.gt = gt; p4.snr_db = 12;
[S4, tau4] = synth_quadpol_column(z, dl_true, deg2rad(35), dt, p4);
o4 = ptt.nodeFabric(S4, z, dt, opts);
th4 = rad2deg(o4.psi0);
e4 = abs(angle(exp(2i * deg2rad(th4 - 35)))) / 2 * 180 / pi;
N4 = o4.nodes;
zn4 = H_ladder(o4, z, tau4, half_wave);
[zerr4, extra4] = H_node_check(N4, zn4, fringe);
ok4m = isfinite(o4.dlam) & o4.zm > 100;
zsc4 = (o4.dlam(ok4m) - dl_t(ok4m)) ./ o4.sigma_dlam(ok4m);
g = o4.gains;
fprintf('\n4. uncalibrated channels, depol 0.5: gains recovered gr %.2f at %.0f deg (truth %.2f at %.0f), gt %.2f at %.0f (truth %.2f at %.0f); firn coherence %.2f\n', ...
  abs(g.gr), rad2deg(angle(g.gr)), abs(gr), rad2deg(angle(gr)), abs(g.gt), rad2deg(angle(g.gt)), abs(gt), rad2deg(angle(gt)), o4.gain_info.coh_firn);
fprintf('   slow axis %.1f (err %.1f), fine_ref %.3f ns, nodes %d, depth err max %.0f m, node resid max %.3f ns\n', ...
  th4, e4, o4.fine_ref, numel(N4.z), zerr4, max(abs(N4.resid_ns)));
fprintf('   found %s\n   true  %s\n', sprintf('%5.0f ', N4.z), sprintf('%5.0f ', zn4));
fprintf('   dlam median |err| %.4f, z rms %.2f, chi2/dof %.2f, abstain below %g m, gamma median %.2f\n', ...
  median(abs(o4.dlam(ok4m) - dl_t(ok4m))), rms(zsc4), o4.chi2_dof, o4.abstain_below_m, median(o4.windows.gamma, 'omitnan'));
o4u = ptt.nodeFabric(S4, z, dt, struct('fc', fc, 'calibrate', false));
ok4u = isfinite(o4u.dlam) & o4u.zm > 100;
fprintf('   for the record, WITHOUT calibration: slow axis %.1f, dlam median |err| %.4f\n', ...
  rad2deg(o4u.psi0), median(abs(o4u.dlam(ok4u) - dl_t(ok4u))));
% the ladder here is the band-centre one. Against it the first two extrema
% sit ~20 m shallow in this column (0.04-0.06 ns, about the node sigma):
% with half the power depolarised the ratio swings little and its extrema
% are pulled up the contrast gradient. 30 m, not the clean column's 15.
ok4 = e4 < 5 && N4.first_ok && zerr4 < 30 && extra4 && median(abs(o4.dlam(ok4m) - dl_t(ok4m))) < 0.005 && ...
  rms(zsc4) < 2.5 && isinf(o4.abstain_below_m);
fprintf('   survives uncalibrated channels and depolarisation:       %s\n', H_tick(ok4));

%% ---- verdict 5: weak fabric, no nodes
dl_weak = 0.003 + 0.004 * (z / 1400);
p5 = p; p5.depol = 0.3;
[S5, tau5] = synth_quadpol_column(z, dl_weak, deg2rad(60), dt, p5);
o5 = ptt.nodeFabric(S5, z, dt, struct('fc', fc, 'psi0', deg2rad(60)));
dl_w = interp1(z, dl_weak, o5.zm + 5, 'linear', 'extrap');
ok5m = isfinite(o5.dlam) & o5.zm > 100;
fprintf('\n5. weak fabric (delay at the bottom %.2f ns): nodes %d, dlam median |err| %.4f, abstain below %g m, inc %d\n', ...
  tau5(end), numel(o5.nodes.z), median(abs(o5.dlam(ok5m) - dl_w(ok5m))), o5.abstain_below_m, o5.n_inc);
if ~isempty(o5.nodes.z)
  fprintf('   spurious node(s): %s\n', sprintf('%.0f m (%+d, %.1f dB, n %g, resid %+.2f ns) ', ...
    [o5.nodes.z.'; o5.nodes.type.'; o5.nodes.prom_db.'; o5.nodes.n.'; o5.nodes.resid_ns.']));
  rs = o5.nodes.smoothed; fprintf('   ratio45 at 100 m steps:'); for zz = 100:100:1400, [~, i] = min(abs(z - zz)); fprintf(' %.1f', rs(i)); end; fprintf('\n');
end
ok5 = (isempty(o5.nodes.z) || ~o5.nodes.first_ok) && median(abs(o5.dlam(ok5m) - dl_w(ok5m))) < 0.004 && isinf(o5.abstain_below_m);
fprintf('   weak fabric from the increments alone:                   %s\n', H_tick(ok5));

%% ---- verdict 6: a trace the cache lacks
S6 = S4;
for ch = {'hh', 'vv', 'hv', 'vh'}, S6.(ch{1})(:, 137) = NaN; end   % the real case: one trace gone in every channel
S6.hv(2000:2100, 251) = NaN;                                        % and a gap in one channel only
o6 = ptt.nodeFabric(S6, z, dt, opts);
ok6m = isfinite(o6.dlam) & isfinite(o4.dlam) & o6.zm > 100;
dth6 = abs(angle(exp(2i * (o6.psi0 - o4.psi0)))) / 2 * 180 / pi;
ddl6 = median(abs(o6.dlam(ok6m) - o4.dlam(ok6m)));
fprintf('\n6. one all-NaN trace and one one-channel gap, traces dropped %d: gains gr %.3f at %.1f deg (was %.3f at %.1f), slow axis moved %.2f deg, dlam median |change| %.5f, nodes %d (was %d)\n', ...
  o6.n_traces_dropped, abs(o6.gains.gr), rad2deg(angle(o6.gains.gr)), abs(g.gr), rad2deg(angle(g.gr)), dth6, ddl6, numel(o6.nodes.z), numel(N4.z));
ok6 = o6.n_traces_dropped == 2 && isfinite(o6.gains.gr) && abs(abs(o6.gains.gr) - abs(g.gr)) < 0.01 && dth6 < 0.5 && ...
  ddl6 < 0.001 && numel(o6.nodes.z) == numel(N4.z) && isinf(o6.abstain_below_m);
fprintf('   a missing trace changes nothing:                         %s\n', H_tick(ok6));

%% ---- verdict 7: the nodes count in the band's frequency
th7 = deg2rad(35);
o7 = struct('fc', fc, 'psi0', th7);
[S7, tau7] = synth_quadpol_column(z, 2 * dl_true, th7, dt, p);
on = ptt.nodeFabric(S7, z, dt, o7);
N7 = on.nodes;
t_at = interp1(z, tau7, N7.z);                          % the true delay where the nodes were found
e_band = rms(t_at - N7.n ./ (2 * N7.fc_hz * 1e-9));
e_lo = rms(t_at - N7.n * fringe / 2);
p7 = p; p7.f_off = 13e6;
S7b = synth_quadpol_column(z, 2 * dl_true, th7, dt, p7);
onb = ptt.nodeFabric(S7b, z, dt, o7);
fprintf('\n7. nine-node column: %d nodes, node frequency %.1f MHz (band centre %.1f, LO %.1f)\n', ...
  numel(N7.z), median(N7.fc_hz) / 1e6, (fc + f_off) / 1e6, fc / 1e6);
fprintf('   true delay at the nodes less n half waves: band %.3f ns rms, LO %.3f ns rms; band moved to +13 MHz -> %.1f MHz\n', ...
  e_band, e_lo, median(onb.nodes.fc_hz) / 1e6);
ok7 = numel(N7.z) >= 8 && abs(median(N7.fc_hz) - (fc + f_off)) < 1.5e6 && e_band < 0.035 && e_lo > 0.05 && ...
  abs(median(onb.nodes.fc_hz) - (fc + 13e6)) < 1.5e6 && isinf(on.abstain_below_m);
fprintf('   the nodes count in the band''s frequency:                 %s\n', H_tick(ok7));

%% ---- verdict 8: a reflection phase is not fabric
p8 = p; p8.refl_phase_deg = 4.5; p8.refl_len_m = 25;
[S8, ~, alpha8] = synth_quadpol_column(z, dl_true, th7, dt, p8);
o8 = ptt.nodeFabric(S8, z, dt, o7);
o8p = ptt.nodeFabric(S8, z, dt, H_merge(o7, poster));
in8 = o8.zm > 100 & o8.zm < 1380;
dl8 = interp1(z, dl_true, o8.zm + 5);                   % cell centres
fine = @(o) rms(o.dlam_raw(in8) - H_pick(movmean(o.dlam_raw, 15), in8));
err8 = rms(o8.dlam(in8) - dl8(in8)); err8p = rms(o8p.dlam(in8) - dl8(in8));
z8 = rms((o8.dlam(in8) - dl8(in8)) ./ o8.sigma_dlam(in8));
% the phase each window saw: the injected one under the window's own taper
wn = o8.windows.win; hw = 0.5 - 0.5 * cos(2 * pi * (0:wn-1).' / (wn - 1));
st = o8.windows.idx - floor(wn / 2);
psi_in = std(arrayfun(@(i) sum(hw .* alpha8(i + (0:wn-1).')) / sum(hw), st)) / (2 * pi * fc * 1e-9);
I8 = o8.intervals;
dlI = arrayfun(@(a, b) mean(dl_true(z >= a & z < b)), I8.z1, I8.z2);
zI = rms((I8.dlam - dlI) ./ I8.sigma, 'omitnan');
fprintf('\n8. reflection phase %.1f deg, %g m: stripes (10 m profile less its 150 m mean) poster %.4f, depth model %.4f\n', ...
  p8.refl_phase_deg, p8.refl_len_m, fine(o8p), fine(o8));
fprintf('   dlam rms error poster %.4f, depth model %.4f (z rms %.2f); resolves %.0f m, %.1f numbers per column\n', ...
  err8p, err8, z8, median(o8.resolution_m(in8), 'omitnan'), o8.n_eff);
fprintf('   window phase found %.4f ns (injected, per window %.4f); residual lag-one correlation %+.2f\n', ...
  o8.depth_model.psi_sd_ns, psi_in, o8.fit.resid_acf1);
fprintf('   prior-free %g m interval means: rms error %.4f, against their stated errors %.2f\n', ...
  I8.z2(1) - I8.z1(1), rms(I8.dlam - dlI, 'omitnan'), zI);
ok8 = fine(o8p) > 0.006 && fine(o8) < 0.002 && err8 < 0.4 * err8p && z8 < 2 && ...
  o8.depth_model.psi_sd_ns > 0.6 * psi_in && o8.depth_model.psi_sd_ns < 1.5 * psi_in && ...
  o8.fit.resid_acf1 < 0 && zI < 2 && rms(I8.dlam - dlI, 'omitnan') < 0.004 && isinf(o8.abstain_below_m);
fprintf('   a reflection phase is not returned as fabric:            %s\n', H_tick(ok8));

%% ---- verdict 9: the answer does not move with the window
wide = struct('win_m', 100, 'step_m', 50);
o9 = ptt.nodeFabric(S8, z, dt, H_merge(o7, wide));
o9p = ptt.nodeFabric(S8, z, dt, H_merge(H_merge(o7, poster), wide));
mv = rms(o8.dlam_raw(in8) - o9.dlam_raw(in8));
mvp = rms(o8p.dlam_raw(in8) - o9p.dlam_raw(in8));
fprintf('\n9. window 50 -> 100 m: the 10 m profile moves %.4f under the depth model (posterior sd %.4f), %.4f under the poster''s\n', ...
  mv, median(o8.sigma_dlam(in8)), mvp);
ok9 = mv < 0.003 && mv < 1.5 * median(o8.sigma_dlam(in8)) && mvp > 2 * mv;
fprintf('   the answer does not move with the window:                %s\n', H_tick(ok9));

%% ---- verdict 10: the answer does not move with the grid
held = H_merge(o7, struct('psi0', o8.psi0, 'gains', o8.gains));
o10 = ptt.nodeFabric(S8, z, dt, held);
cell10 = @(o) accumarray(floor(o.zm(o.zm < 1400) / 20) + 1, o.dlam_raw(o.zm < 1400), [70 1], @mean);
ref10 = cell10(o10);
worst = 0; neff = o10.n_eff;
for g = [20 5 2.5]
  og = ptt.nodeFabric(S8, z, dt, H_merge(held, struct('model_dz', g)));
  worst = max(worst, max(abs(cell10(og) - ref10)));
  neff(end + 1) = og.n_eff; %#ok<SAGROW>
end
fprintf('\n10. grid 10 -> 20, 5, 2.5 m: 20 m means move by at most %.5f (posterior sd %.4f); independent values %s\n', ...
  worst, median(o10.sigma_dlam(in8)), sprintf('%.1f ', neff));
ok10 = worst < 0.1 * median(o10.sigma_dlam(in8)) && max(abs(neff - neff(1))) < 0.5;
fprintf('   the answer does not move with the grid:                  %s\n', H_tick(ok10));

allok = ok1 && ok2 && ok3 && ok4 && ok5 && ok6 && ok7 && ok8 && ok9 && ok10;
fprintf('\n%s  (%.0f s)\n', H_tick(allok), toc(t0));
if ~allok, error('test_node_fabric: FAILED'); end

function [zerr, extra_ok] = H_node_check(N, zn_true, fringe)
% every true node found within its depth error; any extra node must sit
% consistently on the ladder (small residual), e.g. one at the record's edge
zerr = 0; extra_ok = true;
for j = 1:numel(zn_true)
  [dmin, i] = min(abs(N.z - zn_true(j)));
  if isempty(dmin) || dmin > 40 || N.n(i) ~= j, zerr = Inf; else, zerr = max(zerr, dmin); end
end
if numel(N.z) > numel(zn_true)
  extra_ok = all(abs(N.resid_ns(numel(zn_true) + 1:end)) < fringe / 4);
end
end

function zn = H_ladder(o, z, tau, half_wave)
% where the nodes should be: the estimator references every delay to the
% firn window (the channel phase calibration leaves, and the window's own
% delay, go into fine_ref), so node n sits where the true delay below that
% window is n half waves of the band's centre less the firn phase
zw = z(o.windows.idx);
inref = zw >= 30 & zw <= 90 & isfinite(o.windows.fine);
t = tau - interp1(z, tau, mean(zw(inref)));
n = 1:floor((t(end) + o.fine_ref) / half_wave);
zn = arrayfun(@(k) interp1(t, z, k * half_wave - o.fine_ref), n);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function a = H_merge(a, b)
% the fields of b over those of a
fn = fieldnames(b);
for i = 1:numel(fn), a.(fn{i}) = b.(fn{i}); end
end

function v = H_pick(v, mask)
v = v(mask);
end
