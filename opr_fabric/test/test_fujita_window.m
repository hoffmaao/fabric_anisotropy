%TEST_FUJITA_WINDOW The forward model windows like the data path, on any caller grid.
%
% A data path forms its observables at the record's native sampling - the
% 16 moments averaged over a range window W, the co-pol coherence and the
% power anomalies read off those moments - and then keeps one row per W
% so the rows are independent. ptt.fujitaModel used to build its eq.-(7)
% window from however many CALLER rows fell inside win_m, with a floor of
% three, and to take its power anomalies from point amplitudes: handed the
% decimated rows it therefore windowed the coherence over 3W and left the
% powers unwindowed, and model and data disagreed most where the phase
% turns fastest (issue #29). It now evaluates on an internal grid no
% coarser than W/20 and windows there.
%
%   1. AGREEMENT. A one-layer column built on a 0.5 m grid, windowed over
%      W = 20 m with the model's own odd-sample window convention (the
%      one apres.observables uses), decimated to one row per W: the model
%      evaluated at those rows with win_m = win_power_m = W and the data's
%      native step as dz_model reproduces phi, |C| and dP_hh to roundoff
%      (1e-9), because its internal rows then ARE the data's rows and its
%      windows the data's windows. The point-amplitude dP the paper's
%      convention gives (win_power_m = 0) is reported beside it so the
%      size of that mismatch is on record - measured, 5 dB at the nulls.
%   2. CHI-SQUARE AT TRUTH. With Gaussian noise at exactly the sigmas
%      ptt.fabricGLS assigns (the dB variance of an N-look power, the
%      Cramer-Rao phase variance at the observed |C|), the misfit of the
%      TRUE column, standardised by those sigmas, has chi2/dof within
%      0.85-1.15 of one. That is the exit criterion of #29 read at the
%      solution rather than through a solver that may not reach it (#28).
%   3. FINE GRIDS ARE UNTOUCHED. On the 0.5 m grid the refinement factor
%      is 1: out.dz_model is the caller's step, passing that step as
%      dz_model explicitly changes nothing bit for bit, and C and dP_hh
%      equal the pre-change formula - the eq.-(7) window formed by conv2
%      over the CALLER's rows, the anomalies from point amplitudes - to
%      roundoff, computed from a separate evaluation of the stack. The
%      scattering amplitudes of every call are bit-identical to that
%      evaluation's at the same rows, refined or not (verdict 1's
%      41-fold refinement, a 0.1 m grid refined twice), because the
%      internal grid holds the caller's rows as the same doubles; 0.1 is
%      not a binary fraction, so a row rebuilt from a step would not
%      have landed on itself.
%   4. NON-UNIFORM ROWS are interpolated, within the same tolerances as
%      verdict 1, and a bed among them lands on the caller's row that
%      contains it: that row carries the bed's return, every row below
%      it is NaN in every field, the rows clear of the bed's window are
%      unchanged from the bedless call, and a uniform grid with a row at
%      the same depth gives the same bed row - the bed convention is one
%      convention on both grids, evaluated at the caller row's depth.
%   5. ROW ORDER. The same rows given in descending order - on a grid the
%      model refines, on a non-uniform grid, and on one it leaves alone
%      (the 0.5 m grid at dz_model 0.5, refinement factor 1),
%      each with and without a bed - return the ascending call's rows
%      reversed, bit for bit in every per-depth field, with the bed's NaN
%      band still below the bed. The refinement had assumed ascending rows
%      and errored on a descending axis the pre-#29 model had run.
%
% Run: matlab -batch "run('opr_fabric/test/test_fujita_window.m')"
clear;
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
rng(4);

fc = 300e6; W = 20;
DL = 0.23; TH = deg2rad(35); RDB = 0;
dzf = 0.5; zf = (0:dzf:1200).'; Nf = numel(zf);
psi = deg2rad(0:10:170); Np = numel(psi);
lay = struct('top_m', 0, 'dlam', DL, 'theta', TH, 'r_db', RDB, 'gx_db', 0);
opt_fine = struct('fc', fc, 'win_m', 1);          % point values on the fine grid

% --- the data path: moments over W at native sampling, one row per W
fm = ptt.fujitaModel(lay, zf, psi, opt_fine);
nw = 2*floor(W / dzf / 2) + 1; kk = ones(nw, 1) / nw;   % the model's odd window
Phh = conv2(abs(fm.s_hh).^2, kk, 'same');
Pvv = conv2(abs(fm.s_vv).^2, kk, 'same');
Phv = conv2(abs(fm.s_hv).^2, kk, 'same');
Chv = conv2(real(fm.s_hh .* conj(fm.s_vv)), kk, 'same') ...
    + 1i * conv2(imag(fm.s_hh .* conj(fm.s_vv)), kk, 'same');
keep = (nw:nw:Nf - nw).';
z = zf(keep);
A_hh = sqrt(Phh(keep, :)); A_hv = sqrt(Phv(keep, :));
d_hh = 20*log10(A_hh ./ mean(A_hh, 2));
d_hv = 20*log10(max(A_hv, realmin) ./ mean(A_hv, 2));
C = Chv(keep, :) ./ sqrt(Phh(keep, :) .* Pvv(keep, :));
inner = z > 100 & z < 1100;             % clear of the record's edge damping

fails = 0;

% ---- 1. the model on the decimated rows
mo = ptt.fujitaModel(lay, z, psi, struct('fc', fc, 'win_m', W, 'win_power_m', W, 'dz_model', dzf));
e_phi = max(abs(angle(exp(1i * (mo.phi(inner, :) - angle(C(inner, :)))))), [], 'all');
e_cm = max(abs(mo.Cmag(inner, :) - abs(C(inner, :))), [], 'all');
g = isfinite(d_hh(inner, :)) & d_hh(inner, :) > -25;   % where a null's dB is still finite
dd = mo.dP_hh(inner, :) - d_hh(inner, :);
e_hh = max(abs(dd(g)));
mp = ptt.fujitaModel(lay, z, psi, struct('fc', fc, 'win_m', W, 'dz_model', dzf));
dp = mp.dP_hh(inner, :) - d_hh(inner, :);
e_hh_point = max(abs(dp(g)));
ok1 = e_phi < 1e-9 && e_cm < 1e-9 && e_hh < 1e-9 && e_hh_point > 5;
fprintf('1. decimated rows (%d m): phi %.1e rad, |C| %.1e, dP_hh %.1e dB (point-amplitude dP: %.2f dB): %s\n', ...
  W, e_phi, e_cm, e_hh, e_hh_point, H_tick(ok1));
fails = fails + ~ok1;

% ---- 2. chi-square of the true column against noise at the assumed sigmas
NLOOK = 10;
s_db = (10/log(10)) * sqrt(psi_trigamma(NLOOK));
cm = min(max(abs(C), 1e-3), 0.999);
s_ph = sqrt((1 - cm.^2) ./ (2 * NLOOK * cm.^2));
R = 20; chi = nan(R, 1);
for r = 1:R
  o_hh = d_hh + s_db * randn(size(d_hh));
  o_ph = angle(exp(1i * (angle(C) + s_ph .* randn(size(C)))));
  r_hh = (o_hh - mo.dP_hh) / s_db;
  r_ph = angle(exp(1i * (o_ph - mo.phi))) ./ s_ph;
  res = [r_hh(inner, :); r_ph(inner, :)];
  chi(r) = mean(res(:).^2);
end
ok2 = median(chi) > 0.85 && median(chi) < 1.15;
fprintf('2. chi2/dof of the TRUE column at the assumed sigmas: median %.3f over %d draws: %s\n', ...
  median(chi), R, H_tick(ok2));
fails = fails + ~ok2;

% ---- 3. a fine caller grid takes the old code path unchanged
a = ptt.fujitaModel(lay, zf, psi, struct('fc', fc, 'win_m', W));
a_dz = ptt.fujitaModel(lay, zf, psi, struct('fc', fc, 'win_m', W, 'dz_model', dzf));
[e_c3, e_p3] = H_old_path(a, fm, W, dzf);
z1 = (0:0.1:1200).';
a1 = ptt.fujitaModel(lay, z1, psi, struct('fc', fc, 'win_m', W));   % factor 1
f1 = ptt.fujitaModel(lay, z1, psi, opt_fine);                        % factor 2
[e_c1, e_p1] = H_old_path(a1, f1, W, median(diff(z1)));   % the step the old code measured
same_s = H_same_s(a, fm) && H_same_s(a1, f1) && H_same_s(mo, H_rows(fm, keep));
ok3 = a.dz_model == dzf && isequaln(a, a_dz) && same_s ...
  && max([e_c3, e_p3, e_c1, e_p1]) < 1e-9;
fprintf('3. fine grids unrefined: dz_model %.2f, explicit step bit-identical %d, amplitudes bit-identical %d, pre-change formula within %.1e (|C|) %.1e dB: %s\n', ...
  a.dz_model, isequaln(a, a_dz), same_s, max(e_c3, e_c1), max(e_p3, e_p1), H_tick(ok3));
fails = fails + ~ok3;

% ---- 4. non-uniform rows
zn = sort(z + 3 * (rand(size(z)) - 0.5));
mn = ptt.fujitaModel(lay, zn, psi, struct('fc', fc, 'win_m', W, 'win_power_m', W, 'dz_model', dzf));
Cn = interp1(zf, real(Chv ./ sqrt(Phh .* Pvv)), zn, 'linear') ...
   + 1i * interp1(zf, imag(Chv ./ sqrt(Phh .* Pvv)), zn, 'linear');
inn = zn > 100 & zn < 1100;
e_phi4 = max(abs(angle(exp(1i * (mn.phi(inn, :) - angle(Cn(inn, :)))))), [], 'all');
ZB = 1000.3;
mb = ptt.fujitaModel(lay, zn, psi, struct('fc', fc, 'win_m', W, 'win_power_m', W, 'dz_model', dzf, ...
  'bed', struct('z_m', ZB, 'gx_db', 30, 'r_db', 0)));
ib = find(zn >= ZB, 1);
fld = {'s_hh', 's_vv', 's_hv', 'dP_hh', 'dP_hv', 'C', 'phi', 'Cmag', 'P_hh_db'};
fin_above = true; nan_below = true; e_clear = 0;
for f = fld
  fin_above = fin_above && all(isfinite(mb.(f{1})(1:ib, :)), 'all');
  nan_below = nan_below && all(isnan(mb.(f{1})(ib+1:end, :)), 'all');
  e_clear = max(e_clear, max(abs(mb.(f{1})(1:ib-2, :) - mn.(f{1})(1:ib-2, :)), [], 'all'));
end
bed_db = mb.P_hh_db(ib) - max(mb.P_hh_db(1:ib-1));
zu = zn(ib) + (W + dzf) * (-(ib - 1):5).';      % uniform rows through the same bed row
mu = ptt.fujitaModel(lay, zu, psi, struct('fc', fc, 'win_m', W, 'win_power_m', W, 'dz_model', dzf, ...
  'bed', struct('z_m', ZB, 'gx_db', 30, 'r_db', 0)));
iu = find(zu >= ZB, 1);
e_bed = max(max(abs(mb.Cmag(ib, :) - mu.Cmag(iu, :))), ...
  max(abs(angle(exp(1i * (mb.phi(ib, :) - mu.phi(iu, :)))))));
e_bed_db = max(abs(mb.P_hh_db(ib) - mu.P_hh_db(iu)), max(abs(mb.dP_hh(ib, :) - mu.dP_hh(iu, :))));
same_bed = iu == ib && all(isnan(mu.C(iu+1:end, :)), 'all') && e_bed < 0.02 && e_bed_db < 0.1;
okb = fin_above && nan_below && bed_db > 20 && e_clear < 1e-9 && same_bed;
ok4 = e_phi4 < 0.02 && okb;
fprintf('4. non-uniform rows interpolated: phi within %.4f rad; bed row %d at %.1f m: finite to it %d, NaN below %d, +%.1f dB, rows clear of it within %.1e, uniform grid''s bed row within %.4f (|C|, phi) %.3f dB: %s\n', ...
  e_phi4, ib, zn(ib), fin_above, nan_below, bed_db, e_clear, e_bed, e_bed_db, H_tick(ok4));
fails = fails + ~ok4;

% ---- 5. row order
grids = {'refined', z; 'non-uniform', zn; 'fine', zf};
ok5 = true; n5 = 0;
for gi = 1:size(grids, 1)
  za = grids{gi, 2};
  for withbed = [false true]
    o5 = struct('fc', fc, 'win_m', W, 'win_power_m', W, 'dz_model', dzf);
    if withbed, o5.bed = struct('z_m', ZB, 'gx_db', 30, 'r_db', 0); end
    ma = ptt.fujitaModel(lay, za, psi, o5);
    md = ptt.fujitaModel(lay, flipud(za), psi, o5);
    if strcmp(grids{gi, 1}, 'fine')
      ok5 = ok5 && ma.dz_model == dzf && md.dz_model == dzf;
    end
    for f = fld
      ok5 = ok5 && isequaln(md.(f{1}), flipud(ma.(f{1})));
    end
    if withbed
      ok5 = ok5 && all(isnan(md.C(flipud(za) > za(find(za >= ZB, 1)), :)), 'all') ...
        && all(isfinite(md.C(flipud(za) <= za(find(za >= ZB, 1)), :)), 'all');
    end
    n5 = n5 + 1;
  end
end
fprintf('5. descending rows return the ascending call reversed, bit for bit, on %d grid/bed cases: %s\n', ...
  n5, H_tick(ok5));
fails = fails + ~ok5;

fprintf('\n%s (%.1f s)\n', H_tick(fails == 0), toc(t0));
if fails > 0
  error('test_fujita_window:failed', '%d check(s) failed', fails);
end

function y = psi_trigamma(n)
y = psi(1, n);
end

function [e_c, e_p] = H_old_path(m, pt, W, dz)
% the formula before this change, on the caller's rows: the eq.-(7)
% window from however many caller rows fall in W, the anomalies from
% point amplitudes, both from pt's scattering amplitudes at those rows
nw = max(3, 2*floor(W / dz / 2) + 1); kk = ones(nw, 1) / nw;
num = conv2(real(pt.s_hh .* conj(pt.s_vv)), kk, 'same') ...
    + 1i * conv2(imag(pt.s_hh .* conj(pt.s_vv)), kk, 'same');
den = sqrt(conv2(abs(pt.s_hh).^2, kk, 'same') .* conv2(abs(pt.s_vv).^2, kk, 'same'));
C_old = num ./ den;
A = abs(pt.s_hh);
dP_old = 20*log10(A ./ mean(A, 2));
e_c = max(abs(m.C - C_old), [], 'all');
g = dP_old > -25;                       % clear of the nulls' dB amplification
e_p = max(abs(m.dP_hh(g) - dP_old(g)));
end

function ok = H_same_s(a, b)
ok = isequal(a.s_hh, b.s_hh) && isequal(a.s_vv, b.s_vv) && isequal(a.s_hv, b.s_hv);
end

function b = H_rows(a, idx)
b = struct('s_hh', a.s_hh(idx, :), 's_vv', a.s_vv(idx, :), 's_hv', a.s_hv(idx, :));
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
