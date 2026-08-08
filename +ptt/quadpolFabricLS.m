function out = quadpolFabricLS(S, z, opts)
%QUADPOLFABRICLS Fabric from a local least-squares fit of the coherence field.
%
% out = ptt.quadpolFabricLS(S, z, opts)
%
% Fits the measured complex HH-VV coherence C(psi, z), over ALL synthetic
% azimuths at once, to the exact single-column birefringence model in a
% sliding depth window. Three physical parameters per window - the axis
% azimuth theta0, the accumulated two-way phase delta0 at the window
% centre, and its local rate ddelta/dz - plus one nuisance scale gamma for
% the coherence floor. dlam is ddelta/dz converted by the same constant
% ptt.quadpolFabric uses.
%
% WHY THIS EXISTS. ptt.ershadiFabric follows the published chain and takes
% the fabric axis from the cross-polarized power minimum. On this system
% the cross-polarized channels sit on a flat instrument pedestal (cross/co
% ~ -3.6 dB with 0.35 dB of depth ripple, where fabric cross-pol must null
% through every delta = 2*pi*m), so that minimum is antenna-locked: 89.6
% +- 1.9 deg across 1790 blocks spanning all headings. Evaluating the
% phase gradient there, ~25 deg off the true axis, is what produced the
% banded dropouts (odd-pi coherence nulls) and the fringe-periodic
% roughness in the Ridge A section. This estimator never touches the
% cross-polarized power: theta0 comes from the SHAPE of the coherence
% field, which is dominated by the co-polarized channels, and the model
% KNOWS the off-axis nulls, so band samples are fitted rather than gated.
%
% THE MODEL. For a column with principal axes at theta0 and accumulated
% two-way eigenmode phase difference delta(z), the deterministic coherence
% at synthetic azimuth psi is, with mu = cos 2(psi - theta0),
%
%   num = (1-mu^2)/2 + ((1+mu^2)/2) cos delta + i mu sin delta
%   den = 1 - ((1-mu^2)/2) (1 - cos delta)        (= |hh|^2 = |vv|^2)
%   C_model = gamma * num / den
%
% (Rotate diag(e^{i delta}, 1) by (psi - theta0) and form <hh vv*>; den is
% the co-polarized power, whose odd-pi dips are the dropout bands.) Within
% a window delta(z) = delta0 + ddelta*(z - zc). gamma in (0, 1] absorbs
% SNR and volume decorrelation; it is solved in closed form, so the grid
% search is only over (theta0, delta0, ddelta).
%
% AMBIGUITY AND SIGN. (theta0 + 90, -delta0, -ddelta) produces the same
% field, so ddelta >= 0 is enforced and theta0 is unique modulo 180: it is
% the axis along which the (deramped) phase difference grows with depth,
% the same convention ershadiFabric reaches through the polarity of Psi.
% An isotropic window has ddelta = 0 INSIDE the parameter space, so noise
% averages to zero instead of folding to a positive floor the way |Psi|
% does.
%
% TWO-PASS USE. theta0 varies slowly at a divide, so estimate it once per
% frame (opts.theta0 absent), then pass that profile back in for the
% per-block section (opts.theta0 set): the block fits then have two free
% physical parameters and no block-to-block axis jitter.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels hh, vv, hv, vh
%   z     [Nt x 1] depth [m]
%   opts  fc (750e6), psi_step_deg (2), win_short_m (10, coherence
%         multilook), win_fit_m (60, LS window), step_m (15, window
%         centres), dlam_max (0.25), deramped (true),
%         theta0 ([] = estimate; else radians, scalar or [Nw x 1] or a
%         struct('z', 'theta') to interpolate), q_min (0.05, theta0
%         contrast below which theta0 is NaN), dlam_min_theta (0.01,
%         contrast below which theta0 is NaN - a window with no
%         resolvable birefringence has no axes to report, and leakage
%         gives the coherence field azimuth structure of its own that the
%         theta grid would otherwise chase), theta_step_deg (3),
%         weighting ('crb' default: inverse-variance weights from the
%         measured |C|, which suppress the azimuth-selective coherence
%         collapse at the birefringent crossings; 'uniform' to disable)
%
% Output fields (per window centre zw [Nw x 1])
%   zw, theta0, dlam, gamma, resid (weighted rms), q_theta (0..1 contrast
%   of the theta0 cost, 0 = unconstrained), delta0, and dlam_z / theta0_z
%   interpolated back onto z for section wiring. grad_per_dlam echoes the
%   conversion constant.

if nargin < 3, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
psi_step = H_opt(opts, 'psi_step_deg', 2);
win_short = H_opt(opts, 'win_short_m', 10);
win_fit = H_opt(opts, 'win_fit_m', 60);
step_m = H_opt(opts, 'step_m', 15);
dlam_max = H_opt(opts, 'dlam_max', 0.25);
deramped = H_opt(opts, 'deramped', true);
theta0_in = H_opt(opts, 'theta0', []);
q_min = H_opt(opts, 'q_min', 0.05);
dlam_min_theta = H_opt(opts, 'dlam_min_theta', 0.01);
th_step = H_opt(opts, 'theta_step_deg', 3);

z = z(:);
Nt = numel(z);
dz = median(abs(diff(z)));

C = ptt.constants();
n_ice = sqrt(C.eps_bar);
% Same constant as ptt.quadpolFabric, cross-checked there against the
% independent fringe_check.m calibration to 0.7%.
grad_per_dlam = 2 * pi * fc * C.deps / (n_ice * C.c * 1e9);   % rad/m

% --- multilooked coherence over the full sweep. quadpolMoments averages
% the traces; the range term supplies the short window, so A.chhvv is
% already the windowed coherence of Ershadi eq. (7) at every azimuth.
nr = max(3, round(win_short / max(dz, eps)));
M = ptt.quadpolMoments(S, [nr size(S.hh, 2)]);
psi = (0:psi_step:180-psi_step) * pi/180;
A = ptt.quadpolAzimuth(M, psi);
Cm = A.chhvv;
if deramped
  Cm = conj(Cm);
end

% Inverse-variance weights from the measured coherence itself. The
% Cramer-Rao variance of a coherence-phase estimate goes as
% (1 - |C|^2) / (2 N |C|^2), so each (psi, z) sample is weighted by
% |C|^2 / (1 - |C|^2). This suppresses exactly the samples the
% instrument corrupts: at odd-pi crossings the OFF-axis co-polarized
% powers dip into the noise floor, and at even-2pi crossings the
% cross-pol pedestal (decorrelated from the co-pol speckle) enters the
% synthesized T_hh and T_vv with opposite signs; both collapse measured
% |C| azimuth-selectively, and with uniform weights both dragged the fit
% into the gamma ~0.2, resid ~0.6 overshoots the Ridge A validation
% showed at the crossing depths ("regular jumps"). The cap keeps one
% pristine sample from owning a window.
if strcmpi(H_opt(opts, 'weighting', 'crb'), 'crb')
  Wc = abs(Cm).^2 ./ max(1 - abs(Cm).^2, 0.02);
  Wc = min(Wc, 25);
else
  Wc = ones(size(Cm));
end

% --- window centres and in-window sample decimation. ~2 m sampling keeps
% the fit over-determined without dragging 140 correlated samples through
% every grid evaluation.
half = win_fit / 2;
zw = (z(1) + half:step_m:z(end) - half).';
Nw = numel(zw);
jdec = max(1, round(2 / max(dz, eps)));

% theta0 supplied by the caller, in whichever of the three forms
th_fix = nan(Nw, 1);
if ~isempty(theta0_in)
  if isstruct(theta0_in)
    % Interpolate the DOUBLED-ANGLE PHASOR, never an unwrapped angle: an
    % unwrap over gappy, noisy axis samples can slip a branch, which is a
    % silent 90 deg axis error handed to every block below the slip - the
    % 0.25-railed windows at 957 m on frame 007 were exactly that.
    ph = interp1(theta0_in.z(:), exp(2i*theta0_in.theta(:)), zw, ...
      'linear', 'extrap');
    th_fix = 0.5 * angle(ph);
  elseif isscalar(theta0_in)
    th_fix(:) = theta0_in;
  else
    th_fix = theta0_in(:);
    if numel(th_fix) ~= Nw
      error('ptt:quadpolFabricLS:theta0', ...
        'opts.theta0 has %d entries for %d windows', numel(th_fix), Nw);
    end
  end
end

% --- grids. delta0 must cover a full cycle; ddelta >= 0 up to the cap.
th_grid = (0:th_step:180-th_step) * pi/180;
dd_grid = linspace(0, dlam_max * grad_per_dlam, 26);
d0_grid = (0:15:345) * pi/180;

theta0 = nan(Nw, 1); dlam = nan(Nw, 1); gam = nan(Nw, 1);
resid = nan(Nw, 1); q_theta = nan(Nw, 1); delta0 = nan(Nw, 1);

for w = 1:Nw
  jj = find(z >= zw(w) - half & z <= zw(w) + half);
  jj = jj(1:jdec:end);
  Cw = Cm(jj, :).';                 % [Nk x Nj]: azimuth rows, depth cols
  u = (z(jj) - zw(w)).';            % [1 x Nj]
  ok = isfinite(Cw);
  if nnz(ok) < 0.5 * numel(Cw), continue; end
  Cw(~ok) = 0;
  wgt = Wc(jj, :).';
  wgt(~ok) = 0;
  C2 = sum(wgt .* abs(Cw).^2, 'all');
  if C2 <= 0, continue; end

  if isfinite(th_fix(w))
    ths = th_fix(w);
  else
    ths = th_grid;
  end
  best = struct('cost', inf, 'th', NaN, 'd0', NaN, 'dd', NaN);
  cost_th = inf(numel(ths), 1);
  for it = 1:numel(ths)
    mu = cos(2 * (psi(:) - ths(it)));
    Ak = (1 - mu.^2) / 2;
    Bk = (1 + mu.^2) / 2;
    for id = 1:numel(dd_grid)
      cu = cos(dd_grid(id) * u); su = sin(dd_grid(id) * u);
      for i0 = 1:numel(d0_grid)
        c0 = cos(d0_grid(i0)); s0 = sin(d0_grid(i0));
        cd = c0 * cu - s0 * su;
        sd = s0 * cu + c0 * su;
        num = Ak + Bk * cd + 1i * (mu * sd);
        den = max(1 - Ak * (1 - cd), 0.05);
        H = num ./ den;
        G1 = sum(wgt .* real(Cw .* conj(H)), 'all');
        G2 = sum(wgt .* abs(H).^2, 'all');
        cost = C2 - max(G1, 0)^2 / max(G2, realmin);
        if cost < cost_th(it), cost_th(it) = cost; end
        if cost < best.cost
          best = struct('cost', cost, 'th', ths(it), ...
            'd0', d0_grid(i0), 'dd', dd_grid(id));
        end
      end
    end
  end
  if ~isfinite(best.cost), continue; end

  % theta0 contrast: how much worse the fit gets a quarter turn away. An
  % isotropic window is flat here and its theta0 means nothing.
  if numel(ths) > 1
    q_theta(w) = (max(cost_th) - best.cost) / max(C2, realmin);
  end

  % polish from the best grid node. fminsearch is base MATLAB; gamma stays
  % closed-form inside the objective, and ddelta is kept in range by a
  % clamp rather than a penalty cliff. A caller-fixed theta0 stays PINNED:
  % the whole point of the two-pass use is that the block fits cannot
  % wander off the frame axis, so only (delta0, ddelta) are polished then.
  os = optimset('Display', 'off', 'MaxFunEvals', 400, 'MaxIter', 400, ...
    'TolFun', 1e-6, 'TolX', 1e-6);
  fun = @(p) H_cost(p, psi, u, Cw, wgt, C2, dlam_max * grad_per_dlam);
  if isfinite(th_fix(w))
    p2 = fminsearch(@(q) fun([th_fix(w); q(:)]), [best.d0; best.dd], os);
    p = [th_fix(w); p2(:)];
  else
    p = fminsearch(fun, [best.th; best.d0; best.dd], os);
  end
  [cst, g] = fun(p);
  dd = min(max(p(3), 0), dlam_max * grad_per_dlam);

  theta0(w) = mod(p(1), pi);
  delta0(w) = mod(p(2), 2*pi);
  if dd > 0.98 * dlam_max * grad_per_dlam
    % Railed at the cap: the window found no interior optimum, so the
    % value is the bound, not a rate. Abstain rather than report it.
    dlam(w) = NaN;
  else
    dlam(w) = dd / grad_per_dlam;
  end
  gam(w) = g;
  resid(w) = sqrt(max(cst, 0) / C2);
end

if isempty(theta0_in)
  theta0(q_theta < q_min | dlam < dlam_min_theta) = NaN;
end

% interpolate back onto z for the section plumbing; the phasor keeps the
% axis average honest across the modulo-180 wrap
ok = isfinite(theta0);
theta0_z = nan(Nt, 1); dlam_z = nan(Nt, 1);
if nnz(ok) >= 2
  ph = interp1(zw(ok), exp(2i*theta0(ok)), z, 'linear');
  theta0_z = mod(angle(ph)/2, pi);
end
okd = isfinite(dlam);
if nnz(okd) >= 2
  dlam_z = interp1(zw(okd), dlam(okd), z, 'linear');
end

out = struct('zw', zw, 'theta0', theta0, 'dlam', dlam, 'gamma', gam, ...
  'resid', resid, 'q_theta', q_theta, 'delta0', delta0, ...
  'theta0_z', theta0_z, 'dlam_z', dlam_z, ...
  'grad_per_dlam', grad_per_dlam, 'psi', psi);

end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function [cost, gamma] = H_cost(p, psi, u, Cw, wgt, C2, dd_max)
% Weighted misfit of the model against the measured coherence field, with
% the scale gamma eliminated in closed form (clipped to [0, 1.05]: a
% coherence above 1 is not physical, and letting gamma chase one would let
% a noise spike buy a better cost than the data support).
th = p(1);
d0 = p(2);
dd = min(max(p(3), 0), dd_max);
mu = cos(2 * (psi(:) - th));
Ak = (1 - mu.^2) / 2;
Bk = (1 + mu.^2) / 2;
cd = cos(d0 + dd * u);
sd = sin(d0 + dd * u);
num = Ak + Bk * cd + 1i * (mu * sd);
den = max(1 - Ak * (1 - cd), 0.05);
H = num ./ den;
G1 = sum(wgt .* real(Cw .* conj(H)), 'all');
G2 = sum(wgt .* abs(H).^2, 'all');
gamma = min(max(G1, 0) / max(G2, realmin), 1.05);
cost = C2 - 2 * gamma * G1 + gamma^2 * G2;
end
