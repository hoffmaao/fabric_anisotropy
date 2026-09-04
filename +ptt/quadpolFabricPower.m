function out = quadpolFabricPower(S, z, opts)
%QUADPOLFABRICPOWER Fabric from the co-polarized POWER extinction pattern.
%
% out = ptt.quadpolFabricPower(S, z, opts)
%
% Fits the azimuthal pattern of the synthesized co-polarized power
% P(psi, z) in a sliding depth window to the single-column birefringence
% model, with NO use of the HH-VV phase or coherence. It exists for ice
% where that coherence is gone - measured on the Thwaites margin
% (scripts/prototypes/coreg_*_diag.py): the HH-VV pair is depolarised at
% the single-look level, so every estimator built on its phase abstains -
% while the power of each channel still carries the birefringent
% interference of the two eigenmodes, as long as they interfere at all.
%
% THE MODEL. For principal axes at theta0, amplitude reflection
% coefficients r1 (along the axis) and r2 (across it), and accumulated
% two-way differential phase delta(z), the synthesized co-pol amplitude at
% azimuth psi is hh = c^2 r1 e^{i delta} + s^2 r2 with c = cos(psi-theta0),
% s = sin(psi-theta0). With mu = cos 2(psi-theta0), c^2 = (1+mu)/2 and
% s^2 = (1-mu)/2, so the power is
%
%   P = g0 (1+mu)^2/4 + g2 (1-mu)^2/4 + g1 (1-mu^2)/2 cos delta,
%   g0 = r1^2, g2 = r2^2, g1 = r1 r2   (times a common scale)
%
% within a window delta(z) = delta0 + ddelta (z - zc). The model's power
% NODES sit at psi = theta0 +- 45 deg (mu = 0), where P = (g0 + g2)/4 +
% g1/2 cos delta vanishes at delta = pi (2m+1) when r1 = r2: these are the
% extinction nodes of Fujita et al. (2006) and Ershadi et al. (2022), and
% their depth spacing is the birefringence.
%
% A MEASURED LIMITATION (2 Sep 2026, equalised WAIS Divide blocks). The
% node series is one azimuth pair's power against depth, and it
% oscillates for more than one reason: cos delta(z), which is fabric, and
% the reflection ratio r(z) of the layering, which is not. On WAIS
% 20240120_03, where the coherence estimator finds dlam 0.003 with a
% stable axis, the node fit reports |dlam| 0.07-0.14 at visibility
% 0.3-0.45 and explains 62% of the series' variance - against 74-77% on
% Ridge A and WAIS 06 where the fabric is real (dlam 0.06-0.11). The
% r2_min gate therefore does NOT separate the two cases, and a
% power-only |dlam| must not be quoted where the coherence estimator can
% run. What survives at low contrast is the AXIS (mod 90) and the
% visibility; |dlam| from power is for ice where the HH-VV pair is dead
% and the fabric is strong, and it reads ~40% high even there (Ridge A
% 500-800 m, WAIS 06).
%
% WHAT IT CAN AND CANNOT DETERMINE. P depends on delta only through
% cos delta, so the SIGN of ddelta is not observable: dlam is returned as
% a magnitude, and theta0 cannot be assigned to the growing-phase axis
% the way ptt.quadpolFabricLS does. What breaks the 90-deg pair
% {theta0, theta0+90} is the reflection ratio: g0 sits at psi = theta0
% and g2 at theta0 + 90, so with r1 ~= r2 the axis IS unique, and the
% ratio itself comes out (out.r_db). With r1 = r2 the two are degenerate
% and out.theta0 is reported modulo 90 with out.r_db ~ 0 - never
% silently modulo 180.
%
% THE PEDESTAL. The antenna-fixed cross-pol leakage that pedestals the
% coherence (see ptt.quadpolFabricLS) enters the synthesized co-pol power
% as an odd sin 2psi term and an even sin^2 2psi term in the ANTENNA
% frame. Both are linear and are solved with g0, g1, g2 at every node.
% They are collinear with the fabric term when theta0 sits at 0, 45 or 90
% deg to the antennas, exactly as for the coherence; per window they are
% ridged, and the caller can pass a frame-level estimate to anchor them.
%
% All five amplitudes are linear, so at every (theta0, delta0, ddelta)
% grid node the weighted LS reduces to a 5x5 normal system, then the best
% node is polished. Weights are 1/(Pn + w_floor)^2 on the azimuth-
% normalised power: the pattern's information is in its shape, and the
% speckle variance of a power estimate scales with the power itself.
%
% ABSTENTION. A window whose azimuthal modulation depth V = (max-min)/
% (max+min) of the normalised power is below v_min carries no pattern
% (depolarised, isotropic, or noise) and returns NaN; a window whose
% fitted g1^2 disagrees with g0 g2 by more than a factor of q_gg is
% flagged in out.gg_ok (the model's consistency check) but still reported.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels hh, vv, hv, vh; OR a struct
%         with field M ([Nt x 4 x 4] moments), as ptt.quadpolFabricLS
%   z     [Nt x 1] depth [m]
%   opts  fc (750e6), psi_step_deg (2), win_short_m (10, power multilook),
%         win_fit_m (150), step_m (30), dlam_max (0.25),
%         theta_step_deg (3), d0_step_deg (15), n_dd (26),
%         v_min (0.08), w_floor (0.2), tau_rel (0.05, pedestal ridge),
%         q_gg (3), pedestal ('frame' default: two passes, see below;
%         'window': free per window; [a1 a3]: pinned), r_db ([] = fitted,
%         scalar = pinned), r2_min (0.5, variance fraction of the node
%         series the cosine must explain or the window abstains on dlam),
%         r_db_min (0.5, |r_db| below which the 90-deg
%         pair is reported unresolved), ped_az (0; the antenna azimuth in
%         the psi frame, i.e. the heading when S holds geographic moments)
%
% MODE COHERENCE. The interference term's amplitude g1 is fitted freely
% rather than tied to sqrt(g0 g2), and kappa reports the NODE VISIBILITY
% (g1/2) / ((g0+g2)/4 + n0): the swing of the power at the node azimuths
% over its mean, 1 for a fully polarised single column with r1 = r2 and
% falling with depolarised scatter, noise, reflector mismatch between the
% modes and reflection anisotropy. It is the power-only counterpart of
% the HH-VV coherence magnitude (0.3-0.4 on Ridge A 009, where |C| is
% 0.4-0.5), and with n0 (the azimuth-flat unpolarised power) it says from
% power alone how much of a window is fabric signal.
%
% Output (per window centre zw)
%   zw, theta0 (rad, mod pi; see above), dlam (magnitude), delta0,
%   kappa (mode coherence), n0 (unpolarised floor, normalised units),
%   r_db (10 log10 g2/g0: co-pol power across the axis over along it),
%   V (modulation depth), resid (weighted rms of the normalised misfit),
%   q_theta (cost contrast a quarter turn away, on the window's power),
%   gg_ok, ped_coef [Nw x 2], theta0_z / dlam_z on z, grad_per_dlam, psi
%
% See also ptt.quadpolFabricLS, ptt.quadpolAzimuth, ptt.fujitaModel.

if nargin < 3, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
psi_step = H_opt(opts, 'psi_step_deg', 2);
win_short = H_opt(opts, 'win_short_m', 10);
% Longer window than the coherence fit's 60 m: the rate enters here only
% through cos delta, whose sensitivity to ddelta vanishes at every node
% and antinode, so a window must span enough of a cycle to see the
% curvature - 150 m at dlam 0.05 is ~2.2 rad of delta.
win_fit = H_opt(opts, 'win_fit_m', 150);
step_m = H_opt(opts, 'step_m', 30);
dlam_max = H_opt(opts, 'dlam_max', 0.25);
th_step = H_opt(opts, 'theta_step_deg', 3);
d0_step = H_opt(opts, 'd0_step_deg', 15);
n_dd = H_opt(opts, 'n_dd', 26);
v_min = H_opt(opts, 'v_min', 0.08);
w_floor = H_opt(opts, 'w_floor', 0.2);
tau_rel = H_opt(opts, 'tau_rel', 0.05);
q_gg = H_opt(opts, 'q_gg', 3);
% pedestal: 'frame' (default) = two passes, free per window then anchored
% at the frame median with the reflection ratio; 'window' = free per
% window (five amplitudes); [a1 a3] = pinned to a caller value; and
% r_db ([] = fit per window in the free pass, frame median in the anchored
% pass; a scalar pins it)
ped_in = H_opt(opts, 'pedestal', 'frame');
r_db_in = H_opt(opts, 'r_db', []);
% The pedestal is ANTENNA-fixed. When the moments were rotated into a
% geographic frame (ptt.rotateMoments by the heading), the antenna
% azimuth in the psi grid is the heading, and the pedestal shapes are
% sin 2(psi - h) and sin^2 2(psi - h). Ridge A's N-S lines (h ~ 0) hid
% this; WAIS (h ~ 15) and EastGRIP (h ~ 129) do not. Radians.
ped_az = H_opt(opts, 'ped_az', 0);

z = z(:);
Nt = numel(z);
dz = median(abs(diff(z)));
C = ptt.constants();
n_ice = sqrt(C.eps_bar);
grad_per_dlam = 2 * pi * fc * C.deps / (n_ice * C.c * 1e9);   % rad/m

% --- multilooked moments and the synthesized co-pol power
if isfield(S, 'M') && ~isfield(S, 'hh')
  M = S.M;
  nr = max(3, round(win_short / max(dz, eps)));
  kr = ones(nr, 1) / nr;
  for k = 1:4
    for l = 1:4
      M(:, k, l) = conv(M(:, k, l), kr, 'same');
    end
  end
else
  nr = max(3, round(win_short / max(dz, eps)));
  M = ptt.quadpolMoments(S, [nr size(S.hh, 2)]);
end
psi = (0:psi_step:180-psi_step) * pi/180;
A = ptt.quadpolAzimuth(M, psi);
% the two co-pol powers are the same pattern a quarter turn apart
% (T_vv(psi) = T_hh(psi + 90)); using both doubles the looks
P = (A.Phh + A.Pvv(:, H_shift(psi))) / 2;
Pm = mean(P, 2);
Pn = P ./ max(Pm, realmin);                % azimuth-normalised power
Pn(~isfinite(Pn)) = NaN;

half = win_fit / 2;
zw = (z(1) + half:step_m:z(end) - half).';
Nw = numel(zw);
jdec = max(1, round(2 / max(dz, eps)));

th_grid = (0:th_step:180-th_step) * pi/180;
d0_grid = (0:d0_step:360-d0_step) * pi/180;
dd_max = dlam_max * grad_per_dlam;
dd_grid = linspace(0, dd_max, n_dd);
sb = sin(2 * (psi(:) - ped_az)); s2b = sb.^2;   % antenna-frame pedestal shapes

P_ = struct('z', z, 'zw', zw, 'half', half, 'jdec', jdec, 'psi', psi, ...
  'th_grid', th_grid, 'd0_grid', d0_grid, 'dd_grid', dd_grid, 'dd_max', dd_max, ...
  'sb', sb, 's2b', s2b, 'v_min', v_min, 'w_floor', w_floor, 'tau_rel', tau_rel, ...
  'q_gg', q_gg, 'grad_per_dlam', grad_per_dlam, ...
  'r2_min', H_opt(opts, 'r2_min', 0.5));
if ischar(ped_in) && strcmpi(ped_in, 'window')
  R = H_fit_windows(Pn, P_, [], r_db_in);
  R.pedestal = nan(1, 2); R.r_db_frame = NaN;
elseif ischar(ped_in) && strcmpi(ped_in, 'frame')
  % pass A: everything free per window; pass B: pedestal and reflection
  % ratio anchored at their frame medians over the windows with a pattern
  Ra = H_fit_windows(Pn, P_, [], r_db_in);
  okv = isfinite(Ra.V) & Ra.V >= v_min & isfinite(Ra.resid) & isfinite(Ra.theta0);
  if nnz(okv) >= 5
    ped = median(Ra.ped_coef(okv, :), 1);
    if isempty(r_db_in), rdb = median(Ra.r_db(okv)); else, rdb = r_db_in; end
    % EXACT NORMALISATION for pass B. Dividing a row by its azimuth mean
    % (pass A) divides by 3/8 (g0+g2) + g1 cos(delta)/4, which oscillates
    % with delta and bends the normalised pattern away from the linear
    % model - measured as a dlam error of 0.04 on synthetics whose axis
    % was right to a degree. The DC term plus the cos 4(psi-theta0)
    % harmonic is (g0+g2)/2 with NO delta in it, so with the pass-A axis
    % (phasor-interpolated to every row) and the pass-A pedestal removed,
    % that combination normalises every row to the model exactly, and
    % pass B fits the pattern shape with a single scale.
    th_row = H_axis_rows(Ra.theta0, Ra.q_theta, okv, zw, z);
    Pb = P - Pm .* (ped(1) * sb.' + ped(2) * s2b.');
    c4 = 2 * mean(Pb .* cos(4 * (psi - th_row)), 2);
    Nrow = mean(Pb, 2) + c4;
    Pn2 = Pb ./ max(Nrow, realmin);
    Pn2(~isfinite(Pn2) | Nrow <= 0) = NaN;
    % PASS B FITS THE NODE SERIES, not the whole pattern. With the axis
    % known, the power averaged over +-node_half_deg around theta0 +- 45
    % deg is n + a cos(delta(z)) up to a 12% loss of leverage at the band
    % edges - one depth series per row, two linear amplitudes, and the
    % rate enters through cos delta alone. Fitting the full azimuth shape
    % per window instead (tried first) was exact on synthetics and failed
    % on real data wherever the shape departs from the model - depth-
    % varying pedestal cross terms, reflection ratio - by trading depth
    % modulation for shape misfit (dlam 0.001 at 500-800 m on Ridge A
    % 009 against 0.05 from the coherence). The node series is what
    % Fujita and Ershadi read, and it does not care about the rest of the
    % pattern.
    NH = deg2rad(H_opt(opts, 'node_half_deg', 10));
    y = nan(Nt, 1);
    for t = 1:Nt
      if ~isfinite(th_row(t)), continue; end
      dn = abs(angle(exp(2i * (psi - th_row(t) - pi/4)))) / 2;   % distance to theta0+45 mod 90
      dn = min(dn, abs(angle(exp(2i * (psi - th_row(t) + pi/4)))) / 2);
      sel = dn <= NH;
      if any(sel), y(t) = mean(Pn2(t, sel), 'omitnan'); end
    end
    R = H_fit_nodes(y, P_);
    % axis, ratio and pattern diagnostics carry over from pass A
    R.theta0 = Ra.theta0; R.r_db = Ra.r_db; R.V = Ra.V; R.q_theta = Ra.q_theta;
    R.cost_th = Ra.cost_th; R.ped_coef = Ra.ped_coef; R.g = Ra.g; R.gg_ok = Ra.gg_ok;
    R.pedestal = ped; R.r_db_frame = rdb;
    R.Pn2 = Pn2; R.y_node = y;
  else
    R = Ra; R.pedestal = nan(1, 2); R.r_db_frame = NaN;
  end
else
  R = H_fit_windows(Pn, P_, ped_in(:).', r_db_in);
  R.pedestal = ped_in(:).'; R.r_db_frame = r_db_in;
end

% Reflection-ratio degeneracy: where |r_db| is within the noise the pair
% {theta0, theta0+90} is unresolved. Report theta0 modulo 90 there by
% folding to [0, pi/2) so a consumer cannot mistake it for a resolved
% axis; the flag says which windows were folded.
R.theta_mod90 = abs(R.r_db) < H_opt(opts, 'r_db_min', 0.5);
th = R.theta0;
th(R.theta_mod90) = mod(th(R.theta_mod90), pi/2);

theta0_z = nan(Nt, 1); dlam_z = nan(Nt, 1);
okw = isfinite(th);
if nnz(okw) >= 2
  ph = interp1(zw(okw), exp(2i*th(okw)), z, 'linear');
  theta0_z = mod(angle(ph)/2, pi);
end
okd = isfinite(R.dlam);
if nnz(okd) >= 2
  dlam_z = interp1(zw, R.dlam, z, 'linear');
end

out = struct('zw', zw, 'theta0', th, 'theta_mod90', R.theta_mod90, ...
  'dlam', R.dlam, 'delta0', R.delta0, 'r_db', R.r_db, 'V', R.V, ...
  'resid', R.resid, 'q_theta', R.q_theta, 'gg_ok', R.gg_ok, ...
  'g', R.g, 'n0', R.n0, 'kappa', R.kappa, 'ped_coef', R.ped_coef, 'pedestal', R.pedestal, ...
  'node_amp', H_field(R, 'amp', nan(Nw, 1)), 'y_node', H_field(R, 'y_node', nan(Nt, 1)), ...
  'node_r2', H_field(R, 'r2', nan(Nw, 1)), ...
  'r_db_frame', R.r_db_frame, ...
  'theta0_z', theta0_z, 'dlam_z', dlam_z, ...
  'grad_per_dlam', grad_per_dlam, 'psi', psi, 'Pn', Pn);
if isfield(R, 'Pn2'), out.Pn2 = R.Pn2; else, out.Pn2 = Pn; end
end

% -------------------------------------------------------------------------
function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function v = H_field(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end

function th_row = H_axis_rows(th, q, ok, zw, z)
% pass-A axis on every depth row: q-weighted doubled-angle phasor of the
% usable windows, smoothed over 5 windows, interpolated and clamped to the
% window range; rows with nothing usable take the overall phasor mean
qw = max(q, 0); qw(~ok | ~isfinite(qw)) = 0;
ph = qw .* exp(2i * th); ph(~ok) = 0;
ph = conv(ph, ones(5, 1), 'same');
okp = abs(ph) > 0;
if nnz(okp) >= 2
  phz = interp1(zw(okp), ph(okp), min(max(z, min(zw(okp))), max(zw(okp))), 'linear');
else
  phz = repmat(sum(ph), numel(z), 1);
end
th_row = angle(phz) / 2;
end

function idx = H_shift(psi)
% index of psi + 90 deg on the (periodic, uniform) sweep grid
Np = numel(psi);
q = round(Np / 2);
idx = mod((0:Np-1) + q, Np) + 1;
end

function R = H_fit_windows(Pn, P, ped, rdb)
%H_FIT_WINDOWS One pass of per-window fits. ped = [] frees the two pedestal
% amplitudes, [a1 a3] pins them; rdb = [] frees g0, g2, g1 (three
% amplitudes), a scalar ties them to one scale g with g2 = g R^2, g1 = g R.
Nw = numel(P.zw);
Nth = numel(P.th_grid);
R = struct('theta0', nan(Nw,1), 'dlam', nan(Nw,1), 'delta0', nan(Nw,1), ...
  'r_db', nan(Nw,1), 'V', nan(Nw,1), 'resid', nan(Nw,1), ...
  'q_theta', nan(Nw,1), 'gg_ok', false(Nw,1), 'ped_coef', nan(Nw,2), ...
  'g', nan(Nw,3), 'n0', nan(Nw,1), 'kappa', nan(Nw,1), 'cost_th', nan(Nw, Nth));
os = optimset('Display', 'off', 'MaxFunEvals', 400, 'MaxIter', 400, ...
  'TolFun', 1e-7, 'TolX', 1e-6);
Rlin = [];
if ~isempty(rdb), Rlin = 10^(rdb / 20); end
for w = 1:Nw
  jj = find(P.z >= P.zw(w) - P.half & P.z <= P.zw(w) + P.half);
  jj = jj(1:P.jdec:end);
  Y = Pn(jj, :).';                          % [Npsi x Nj]
  ok = isfinite(Y);
  if nnz(ok) < 0.5 * numel(Y), continue; end
  u = (P.z(jj) - P.zw(w)).';
  ybar = mean(Y, 2, 'omitnan');
  V = (max(ybar) - min(ybar)) / max(max(ybar) + min(ybar), realmin);
  R.V(w) = V;
  if ~(V >= P.v_min), continue; end
  Wt = 1 ./ (Y + P.w_floor).^2;
  Wt(~ok) = 0; Y(~ok) = 0;
  Y2 = sum(Wt .* Y.^2, 'all');
  if Y2 <= 0, continue; end
  tau = P.tau_rel * sum(Wt, 'all') * mean(P.s2b);
  K = struct('sb', P.sb, 's2b', P.s2b, 'tau', tau, 'ped', ped, 'R', Rlin, ...
    'dd_max', P.dd_max);

  best = struct('cost', inf, 'th', NaN, 'd0', NaN, 'dd', NaN);
  cost_th = inf(Nth, 1);
  for it = 1:Nth
    mu = cos(2 * (P.psi(:) - P.th_grid(it)));
    F = H_shapes(mu);
    for id = 1:numel(P.dd_grid)
      cu = cos(P.dd_grid(id) * u); su = sin(P.dd_grid(id) * u);
      for i0 = 1:numel(P.d0_grid)
        cd = cos(P.d0_grid(i0)) * cu - sin(P.d0_grid(i0)) * su;
        cst = H_solve(Y, Wt, Y2, F, cd, K);
        if cst < cost_th(it), cost_th(it) = cst; end
        if cst < best.cost
          best = struct('cost', cst, 'th', P.th_grid(it), 'd0', P.d0_grid(i0), ...
            'dd', P.dd_grid(id));
        end
      end
    end
  end
  if ~isfinite(best.cost), continue; end
  R.q_theta(w) = (max(cost_th) - best.cost) / Y2;
  R.cost_th(w, :) = cost_th(:).';

  fun = @(p) H_cost(p, P.psi, u, Y, Wt, Y2, K);
  p = fminsearch(fun, [best.th; best.d0; best.dd], os);
  [cst, x] = fun(p);
  dd = min(max(p(3), 0), P.dd_max);
  R.theta0(w) = mod(p(1), pi);
  R.delta0(w) = mod(p(2), 2*pi);
  if dd > 0.98 * P.dd_max
    R.dlam(w) = NaN;                        % railed at the cap: not a rate
  else
    R.dlam(w) = dd / P.grad_per_dlam;
  end
  R.g(w, :) = x(1:3).';
  R.ped_coef(w, :) = x(4:5).';
  R.n0(w) = x(6);
  % node visibility: the swing of the node series, g1/2, over its mean,
  % (g0+g2)/4 + n0. Bounded in [0, 1] for physical amplitudes and immune
  % to the DC trade between the pattern scale and the unpolarised floor
  % (the ratio g1/sqrt(g0 g2) is not: it read 2 on Ridge A).
  R.kappa(w) = (x(2) / 2) / max((x(1) + x(3)) / 4 + max(x(6), 0), realmin);
  R.resid(w) = sqrt(max(cst, 0) / Y2);
  g0 = max(x(1), realmin); g2 = max(x(3), realmin); g1 = x(2);
  R.r_db(w) = 10 * log10(g2 / g0);
  R.gg_ok(w) = g1 > 0 && g1^2 / (g0 * g2) < P.q_gg && g1^2 / (g0 * g2) > 1 / P.q_gg;
end
end

function R = H_fit_nodes(y, P)
%H_FIT_NODES Per-window fit of the node series y(z) = n + a cos(d0 + dd u).
% Grid over (d0, dd) with the two amplitudes in closed form (a >= 0), then
% a polish. kappa = a / n is the node visibility; resid is the rms misfit
% relative to n.
Nw = numel(P.zw);
R = struct('dlam', nan(Nw,1), 'delta0', nan(Nw,1), 'kappa', nan(Nw,1), ...
  'n0', nan(Nw,1), 'resid', nan(Nw,1), 'amp', nan(Nw,1), 'r2', nan(Nw,1));
os = optimset('Display', 'off', 'MaxFunEvals', 300, 'MaxIter', 300, ...
  'TolFun', 1e-8, 'TolX', 1e-6);
for w = 1:Nw
  jj = find(P.z >= P.zw(w) - P.half & P.z <= P.zw(w) + P.half);
  jj = jj(1:P.jdec:end);
  yy = y(jj); ok = isfinite(yy);
  if nnz(ok) < 0.5 * numel(yy) || nnz(ok) < 8, continue; end
  yy = yy(ok); u = P.z(jj(ok)) - P.zw(w);
  best = struct('cost', inf, 'd0', NaN, 'dd', NaN, 'b', [NaN NaN]);
  for id = 1:numel(P.dd_grid)
    cu = cos(P.dd_grid(id) * u); su = sin(P.dd_grid(id) * u);
    for i0 = 1:numel(P.d0_grid)
      cd = cos(P.d0_grid(i0)) * cu - sin(P.d0_grid(i0)) * su;
      [cst, b] = H_lin2(yy, cd);
      if cst < best.cost
        best = struct('cost', cst, 'd0', P.d0_grid(i0), 'dd', P.dd_grid(id), 'b', b);
      end
    end
  end
  if ~isfinite(best.cost), continue; end
  fun = @(p) H_lin2(yy, cos(p(1) + min(max(p(2), 0), P.dd_max) * u));
  p = fminsearch(fun, [best.d0; best.dd], os);
  [cst, b] = fun(p);
  dd = min(max(p(2), 0), P.dd_max);
  % NULL TEST. A cosine with free phase and rate fits SOME oscillation into
  % any series, so the rate is reported only where it explains a set
  % fraction of the variance the constant-only model leaves: measured on
  % WAIS 03, a near-isotropic frame (coherence dlam 0.003), the ungated
  % fit reported dlam 0.07-0.14 at visibility 0.3-0.45. r2_min is the
  % fraction of the residual variance about the mean removed by the
  % cosine.
  c0 = sum((yy - mean(yy)).^2);
  R.r2(w) = 1 - cst / max(c0, realmin);
  if R.r2(w) < P.r2_min
    R.kappa(w) = b(2) / max(b(1), realmin);
    R.n0(w) = b(1); R.amp(w) = b(2);
    R.resid(w) = sqrt(cst / numel(yy)) / max(b(1), realmin);
    continue;                       % dlam and delta0 stay NaN: abstained
  end
  R.delta0(w) = mod(p(1), 2*pi);
  if dd > 0.98 * P.dd_max
    R.dlam(w) = NaN;
  else
    R.dlam(w) = dd / P.grad_per_dlam;
  end
  R.n0(w) = b(1); R.amp(w) = b(2);
  R.kappa(w) = b(2) / max(b(1), realmin);
  R.resid(w) = sqrt(cst / numel(yy)) / max(b(1), realmin);
end
end

function [cost, b] = H_lin2(yy, cd)
% y = b1 + b2 cd with b2 >= 0
X = [ones(size(cd)), cd(:)];
if std(cd) < 1e-9
  % dd = 0 (or a window too short to see it): the cosine column is a
  % constant and the fit is the mean alone
  b = [mean(yy); 0];
else
  b = X \ yy;
  if b(2) < 0
    b = [mean(yy); 0];
  end
end
cost = sum((yy - X * b).^2);
end

function F = H_shapes(mu)
% the three fabric azimuth shapes for g0, g1, g2: columns [f0 f1 f2]
F = [(1 + mu).^2 / 4, (1 - mu.^2) / 2, (1 - mu).^2 / 4];
end

function [cost, x] = H_solve(Y, Wt, Y2, F, cd, K)
% Weighted LS for the linear amplitudes given the azimuth shapes F =
% [f0 f1 f2] [Npsi x 3], the depth factor cd [1 x Nj] on the f1 term and
% the two antenna-frame pedestal shapes. Amplitude vector x = [g0 g1 g2
% a1 a3]. Three regimes: all five free; pedestal pinned (x(4:5) given);
% reflection ratio tied (g1 = g R, g2 = g R^2 -> one fabric amplitude g).
% Every normal system carries a tiny ridge on its whole diagonal, because
% at theta0 = 0 or 90 deg to the antennas the f1 and sin^2 2psi shapes are
% exactly collinear and the pedestal ridge alone leaves it singular.
wrow = sum(Wt, 2);
wcd = Wt * cd(:);
wcd2 = Wt * (cd(:).^2);
yrow = sum(Wt .* Y, 2);
ycd = sum(Wt .* Y .* cd, 2);
f0 = F(:,1); f1 = F(:,2); f2 = F(:,3);
% UNPOLARISED FLOOR. Additive noise and depolarised scatter add power
% that is flat over azimuth, so the nodes never reach the model's zero.
% Without a term for it the fit - whose node depth is tied to the
% reflection ratio - compromises by flattening its DEPTH modulation, and
% dlam comes out ~0.4x on noisy synthetics (exact without noise). A
% constant column n0 >= 0 per window carries that power; it is
% collinear with f0 + f1 + f2 = 1 only while cd is constant across the
% window, which the ridge covers.
one = ones(size(f0));
if isempty(K.R)
  % basis: f0, f1*cd, f2, sb, s2b, 1 -> x = [g0 g1 g2 a1 a3 n0]
  Bc = [f0, f2, K.sb, K.s2b, one];           % depth-constant columns
  G = zeros(6); r = zeros(6, 1);
  ic = [1 3 4 5 6];
  G(ic, ic) = Bc.' * (Bc .* wrow);
  G(2, 2) = f1.' * (f1 .* wcd2);
  G(2, ic) = (f1 .* wcd).' * Bc;
  G(ic, 2) = G(2, ic).';
  r(ic) = Bc.' * yrow;
  r(2) = f1.' * ycd;
  ip = [4 5]; ineg = [1 2 3 6];
else
  % tied reflection ratio, FREE interference amplitude. The two eigenmode
  % echoes need not interfere with full visibility: partially coherent
  % modes (different reflector sets, depolarised scatter, noise) give a
  % cos delta swing of kappa r1 r2 with kappa < 1 - measured ~0.3 on Ridge
  % A 009, where the HH-VV coherence is ~0.4-0.5 for the same reason.
  % Tying g1 to sqrt(g0 g2) forces a full swing and the fit pays for it
  % with a flattened depth modulation (dlam 3x low at 500-1150 m). So:
  % basis h0 = f0 + R^2 f2 (scale g), f1 cd (scale g1 = kappa g R), sb,
  % s2b, 1 -> x = [g g1 a1 a3 n0], expanded to [g0 g1 g2 a1 a3 n0]
  Bc = [f0 + K.R^2 * f2, K.sb, K.s2b, one];
  G = zeros(5); r = zeros(5, 1);
  ic = [1 3 4 5];
  G(ic, ic) = Bc.' * (Bc .* wrow);
  G(2, 2) = f1.' * (f1 .* wcd2);
  G(2, ic) = (f1 .* wcd).' * Bc;
  G(ic, 2) = G(2, ic).';
  r(ic) = Bc.' * yrow;
  r(2) = f1.' * ycd;
  ip = [3 4]; ineg = [1 2 5];
end
n = numel(r);
free = true(n, 1);
G = G + (1e-9 * trace(G) / n + realmin) * eye(n);
G(ip(1), ip(1)) = G(ip(1), ip(1)) + K.tau;
G(ip(2), ip(2)) = G(ip(2), ip(2)) + K.tau;
x = zeros(n, 1);
if ~isempty(K.ped)
  x(ip) = K.ped(:);
  free(ip) = false;
end
x(free) = G(free, free) \ (r(free) - G(free, ~free) * x(~free));
% physical amplitudes are non-negative: clamp and refit the rest
neg = false(n, 1); neg(ineg) = free(ineg) & x(ineg) < 0;
if any(neg)
  x(neg) = 0; free(neg) = false;
  x(free) = G(free, free) \ (r(free) - G(free, ~free) * x(~free));
end
if ~all(isfinite(x))
  cost = Y2; x = nan(6, 1); return;
end
cost = Y2 - 2 * (x.' * r) + x.' * G * x;
if ~isempty(K.R)
  x = [x(1); x(2); K.R^2 * x(1); x(3); x(4); x(5)];
end
end

function [cost, x] = H_cost(p, psi, u, Y, Wt, Y2, K)
th = p(1); d0 = p(2); dd = min(max(p(3), 0), K.dd_max);
mu = cos(2 * (psi(:) - th));
cd = cos(d0 + dd * u);
[cost, x] = H_solve(Y, Wt, Y2, H_shapes(mu), cd, K);
end

