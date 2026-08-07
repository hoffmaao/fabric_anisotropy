function out = quadpolFabric(A, z, opts)
%QUADPOLFABRIC Horizontal fabric orientation and contrast from a quad-pol sweep.
%
% out = ptt.quadpolFabric(A, z, opts)
%
% Inverts the synthetic azimuth sweep from ptt.quadpolAzimuth for
%
%   out.theta  [Nt x 1]  azimuth of a horizontal principal axis [rad],
%                        modulo pi/2 (see AMBIGUITY below)
%   out.dlam   [Nt x 1]  lam_max - lam_min from the coherence phase gradient
%   out.dlam_node [Nt x 1] the same from cross-polarized node spacing
%   out.delta  [Nt x 1]  unwrapped differential phase at the aligned azimuth
%   out.aniso  [Nt x 1]  0..1 depth of the cross-pol azimuthal modulation
%
% ORIENTATION. The cross-polarized power of a birefringent column goes as
%   Pxc(psi) ~ sin^2(2(psi - theta)) = 1/2 - 1/2 cos(4(psi - theta)),
% i.e. a pure 4-psi harmonic, so theta is obtained by projecting the sweep
% onto cos 4psi and sin 4psi rather than by searching for a minimum. The
% projection uses every azimuth, so it is not thrown by the single deepest
% sample the way argmin is, and its amplitude relative to the mean gives
% `aniso`, a direct measure of how anisotropic that depth actually is - a
% near-isotropic layer returns a flat sweep and a meaningless theta, which
% argmin would report just as confidently as a real one.
%
% AMBIGUITY. The 4-psi harmonic has period pi/2, so this fixes theta only
% modulo 90 degrees: it says where the principal axes are, not which of the
% two carries lam_max. That is broken by the SIGN of the coherence phase
% gradient, which is positive along one axis and negative along the other.
% The convention here returns theta as the axis along which delta INCREASES
% with depth, and flags the choice in out.branch_flipped.
%
% CONTRAST. Two independent estimates, deliberately:
%   dlam       from d(delta)/dz, the unwrapped coherence phase at the
%              aligned azimuth - uses every depth sample, but inherits any
%              unwrapping slip
%   dlam_node  from the depth spacing of the cross-polarized nulls, where
%              delta passes 2*pi*m - immune to unwrapping, but coarse,
%              since it only reports once per fringe
% They fail in unrelated ways, so their agreement is a reliability test
% neither provides alone. Both are the TRUE contrast, not the projected
% one: delta is the phase difference between the medium's own eigenmodes,
% which the antenna azimuth cannot change.
%
% opts fields (all optional)
%   fc        centre frequency [Hz], default 750e6
%   win_m     depth window for the azimuth projection [m], default 50
%   grad_win_m window for the d(delta)/dz slope [m], default 200
%   aniso_min modulation depth below which theta is set NaN, default 0.05

if nargin < 3, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
win_m = H_opt(opts, 'win_m', 50);
grad_win_m = H_opt(opts, 'grad_win_m', 200);
aniso_min = H_opt(opts, 'aniso_min', 0.05);

z = z(:);
Nt = numel(z);
if size(A.Pxc, 1) ~= Nt
  error('ptt:quadpolFabric:size', ...
    'z has %d entries but the sweep has %d rows', Nt, size(A.Pxc,1));
end
psi = A.psi(:).';
dz = median(abs(diff(z)));
nw = max(3, 2*floor(win_m / max(dz, eps) / 2) + 1);

C = ptt.constants();
n_ice = sqrt(C.eps_bar);
% d(delta)/dz per unit dlam. Two-way delay difference accumulates at
% d(dtau)/dz = 2*dn/c with dn = deps*dlam/(2n), and delta = 2*pi*fc*dtau.
% Cross-checked against the independent fringe_check.m calibration
% (4.021 cycles/us per unit dlam at 750 MHz) in egrip_zeising.m: 6.35e-11
% against 6.39e-11 s/m, agreeing to 0.7%.
grad_per_dlam = 2 * pi * fc * 2 * (C.deps / (2*n_ice)) / (C.c * 1e9);

% --- orientation: project the cross-pol sweep onto the 4-psi harmonic
Pxc = movmean(A.Pxc, nw, 1, 'omitnan');
m0 = mean(Pxc, 2);
c4 = cos(4*psi);
s4 = sin(4*psi);
% Pxc = m0 - K/2 * cos(4(psi-theta)) => a1 = -K/2 cos4theta, b1 = -K/2 sin4theta
a1 = 2 * mean(Pxc .* c4, 2);
b1 = 2 * mean(Pxc .* s4, 2);
aniso = hypot(a1, b1) ./ max(m0, realmin);

% Smooth the PHASOR, not the angle. theta is only defined modulo pi/2, so
% a fabric axis that happens to sit near that boundary - which is exactly
% what an antenna pair aligned with the fabric produces - makes the
% per-depth angle jitter between ~0 and ~90 deg under noise. The coherence
% phase is +delta along one of those and -delta along the other, so the
% jitter alternates the sign of the phase profile and destroys the unwrap
% below it. In the synthetic test that read dlam = 0.106 against a truth of
% 0.300, while the raw column at the correct azimuth gave 0.300 exactly.
%
% -(a1 + i b1) is K/2 * exp(i 4 theta), so smoothing it and taking the
% argument averages the ORIENTATION correctly across the wrap, and weights
% each depth by its own modulation amplitude for free. Fabric orientation
% varies slowly with depth, so this is a physical smoothing, not cosmetic.
P4 = -(a1 + 1i*b1);
P4(~isfinite(P4)) = 0;
P4(aniso < aniso_min) = 0;      % flat sweeps carry no orientation
ns = max(3, 2*floor(grad_win_m / max(dz, eps) / 2) + 1);
P4s = movmean(P4, ns, 1, 'omitnan');
if ~any(P4s)
  out = H_empty(Nt, grad_per_dlam);
  return;
end
% UNWRAP the quadruple angle down the column before dividing by four.
%
% theta is only defined modulo pi/2, so every folded representation has a
% branch cut somewhere and an axis sitting on that cut jitters across it
% under noise. The coherence phase is +delta on one side and -delta on the
% other, so the jitter alternates the sign of the phase profile and
% destroys the unwrap below it. Both fixes tried before this one merely
% MOVED the cut: folding to [0, pi/2) failed for antennas aligned with the
% fabric (dlam 0.134 against a truth of 0.300), and taking the raw
% argument moved the failure to axes at 45 deg (0.152). Unwrapping removes
% the cut instead of relocating it, and is valid because P4s has already
% been smoothed over grad_win_m and fabric orientation varies slowly.
th_use = 0.25 * unwrap(angle(P4s));
theta = mod(th_use, pi/2);
theta(~isfinite(aniso) | aniso < aniso_min) = NaN;

% --- differential phase at the aligned azimuth, by interpolating the
% sweep rather than snapping to the nearest sampled psi
Chhvv = H_interp_psi(A.Chhvv, psi, th_use);
delta = H_unwrap_finite(angle(Chhvv));

% The 90 deg branch: along one principal axis delta increases with depth,
% along the other it decreases. Take the increasing one, so dlam is
% reported positive for the axis carrying lam_max.
ok = isfinite(delta);
branch_flipped = false;
if nnz(ok) > 8
  p = polyfit(z(ok), delta(ok), 1);
  if p(1) < 0
    branch_flipped = true;
    % th_use stays continuous here for the same reason it is not folded
    % above; only the reported angle is wrapped.
    th_use = th_use + pi/2;
    theta = mod(th_use, pi/2);
    theta(~isfinite(aniso) | aniso < aniso_min) = NaN;
    Chhvv = H_interp_psi(A.Chhvv, psi, th_use);
    delta = H_unwrap_finite(angle(Chhvv));
  end
end

% --- contrast from the phase gradient
ng = max(3, 2*floor(grad_win_m / max(dz, eps) / 2) + 1);
half = floor(ng/2);
grad = nan(Nt, 1);
for i = 1:Nt
  lo = max(1, i-half); hi = min(Nt, i+half);
  zz = z(lo:hi); yy = delta(lo:hi);
  m = isfinite(yy);
  if nnz(m) < 5, continue; end
  zz = zz(m); yy = yy(m);
  if (max(zz) - min(zz)) < 0.5*grad_win_m, continue; end
  pf = polyfit(zz, yy, 1);
  grad(i) = pf(1);
end
dlam = grad / grad_per_dlam;

% --- contrast from cross-polarized node spacing, at 45 deg to the axes
% where the modulation is deepest and the nulls are cleanest
Pnode = H_interp_psi(A.Pxc, psi, mod(th_use + pi/4, pi));
dlam_node = H_node_rate(Pnode, z, grad_per_dlam, grad_win_m);

out = struct('theta', theta, 'dlam', dlam, 'dlam_node', dlam_node, ...
  'delta', delta, 'aniso', aniso, 'grad', grad, ...
  'branch_flipped', branch_flipped, 'grad_per_dlam', grad_per_dlam);

end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function y = H_interp_psi(X, psi, th)
% Value of X(t, :) at azimuth th(t), by linear interpolation on the
% periodic psi grid. Sampling the nearest column instead would quantise
% theta to the sweep spacing, which at 2 deg is a tenth of the difference
% the EastGRIP branch test has to resolve.
Nt = size(X, 1);
y = nan(Nt, 1);
if ~isreal(X), y = complex(y); end
P = psi(end) - psi(1) + (psi(2) - psi(1));   % full period of the grid
for i = 1:Nt
  if ~isfinite(th(i)), continue; end
  u = mod(th(i) - psi(1), P) / (psi(2) - psi(1));
  k = floor(u);
  f = u - k;
  i0 = mod(k, numel(psi)) + 1;
  i1 = mod(k+1, numel(psi)) + 1;
  y(i) = (1-f) * X(i, i0) + f * X(i, i1);
end
end

function dl = H_node_rate(P, z, grad_per_dlam, win_m)
% Local fringe rate of the cross-polarized power, from its SPECTRUM in a
% sliding depth window rather than from the spacing of individual minima.
%
% Pxc ~ sin^2(delta/2) = (1 - cos delta)/2, so it oscillates at the local
% fringe rate and the rate is the dominant frequency of that oscillation.
% Hunting for minima instead was the first implementation and it failed at
% every azimuth of the synthetic test, reading 1.2-1.8 against a truth of
% 0.30: at 0.5 m sampling and realistic speckle there are many spurious
% local minima between the true ones, every one of which shortens the
% apparent spacing and inflates the rate. A spectral peak is set by the
% periodicity of the whole window and is unmoved by them.
%
% This estimator never unwraps a phase, so it fails in a different way
% from the coherence-gradient estimate and their agreement is a real test.
Nt = numel(z);
dl = nan(Nt, 1);
p = P(:);
dz = median(abs(diff(z)));
nw = max(16, 2*floor(win_m / max(dz, eps) / 2));
if nnz(isfinite(p)) < nw, return; end
half = floor(nw/2);
w = hann(nw);
% Frequencies resolvable in this window, skipping DC and the first bin
% (a linear trend across the window, not a fringe).
f = (0:nw-1)' / (nw * dz);
kmax = floor(nw/2);
for i = 1:Nt
  lo = i - half; hi = lo + nw - 1;
  if lo < 1 || hi > Nt, continue; end
  seg = p(lo:hi);
  if any(~isfinite(seg)), continue; end
  seg = seg - mean(seg);
  if ~any(seg), continue; end
  Fv = abs(fft(seg .* w));
  Fv(1:2) = 0;                       % DC and the window-scale trend
  [~, k] = max(Fv(1:kmax));
  dl(i) = 2*pi*f(k) / grad_per_dlam;
end
end

function y = H_unwrap_finite(x)
% unwrap() propagates NaN through everything after the first gap, which is
% fatal here because the sweep goes flat once per fringe by construction.
% Unwrap over the finite samples only and put the gaps back afterwards.
y = nan(size(x));
m = isfinite(x);
if nnz(m) < 2, return; end
y(m) = unwrap(x(m));
end

function out = H_empty(Nt, g)
z = nan(Nt, 1);
out = struct('theta', z, 'dlam', z, 'dlam_node', z, 'delta', z, ...
  'aniso', z, 'grad', z, 'branch_flipped', false, 'grad_per_dlam', g);
end
