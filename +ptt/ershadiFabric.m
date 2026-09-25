function out = ershadiFabric(S, z, opts)
%ERSHADIFABRIC Ershadi et al. (2022) polarimetric fabric algorithm.
%
% out = ptt.ershadiFabric(S, z, opts)
%
% Ershadi, M. R., Drews, R., Martin, C., Eisen, O., Ritz, C., Corr, H.,
% Christmann, J., Zeising, O., Humbert, A. and Mulvaney, R. (2022),
% "Polarimetric radar reveals the spatial distribution of ice fabric at
% domes and divides in East Antarctica", The Cryosphere 16, 1719-1739.
%
% A faithful implementation of the published chain, kept SEPARATE from
% ptt.quadpolFabric so the two can be run on the same data and compared.
% Where this file departs from that one it is because the paper does, and
% each departure is marked (E1)-(E5) below.
%
%   (E1) SYNTHESIS SENSE. Their eq. (11) is
%          S_N(theta +- gamma) = R(theta +- gamma) S_N(theta) R'(...)
%        i.e. R S R'. ptt.rotatePolarization applies R' S R, the inverse
%        rotation, so the two recover orientations of opposite sign
%        (theta -> -theta, which for a modulo-180 axis is 180 - theta).
%        This follows the paper.
%
%   (E2) DERAMPING. "The deramping corresponds to a complex conjugation of
%        C_HHVV", so the paper uses eq. (7) for models and its CONJUGATE
%        for radar data. That flips the sign of phi_HHVV, hence of Psi and
%        of the branch that decides which axis carries lam_max. Controlled
%        by opts.deramped, default true for real data.
%
%   (E3) UNWRAP-FREE GRADIENT. They follow Jordan et al. (2019) and take a
%        1D convolutional derivative of the real and imaginary parts of the
%        COMPLEX coherence, never of the phase:
%          dphi/dz = Im(conj(C) dC/dz) / |C|^2
%        This is not a refinement, it removes a failure mode: differencing
%        an unwrapped phase inherits every cycle slip, and unwrapping is
%        what produced three separate branch-cut bugs in quadpolFabric.
%
%   (E4) ORIENTATION FROM THE POWER ANOMALY. Their eq. (12) works in dB
%        against the azimuthal mean,
%          dP_xx(theta,z) = 20 log10( |s_xx| / mean_b |s_xx(theta_b)| )
%        and orientation is the MINIMUM of dP_HV (cross-polarization
%        extinction), disambiguated by the polarity of phi_HHVV.
%
%   (E5) QUALITY GATE. "We restrict the analysis to the top 2000 m, where
%        typically |C_HHVV| > 0.4."
%
% CONSTANT. Their eq. (9) as printed reads Psi = 2 c eps' /(4 pi fc
% Deps') dphi/dz. Deriving it: the two-way phase difference accumulates at
% dphi/dz = 2 (2 pi fc/c) Dn with Dn = Deps' dlam / (2 sqrt(eps')), giving
%   dlam = c sqrt(eps') / (2 pi fc Deps') dphi/dz,
% i.e. sqrt(eps'), not eps'. The two differ by a factor sqrt(3.15) = 1.77.
% The square-root form is used here because it is what the physics gives,
% because it agrees to 0.7% with this repo's independent calibration of
% ptt.twttDifference (fringe_check.m), and because a 1.77x error would
% have been obvious against the ice cores the paper validates on - so the
% printed form is most likely a typesetting loss of the radical.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels hh, vv, hv (vh optional).
%         A trace may be NaN below its own bed (ptt.maskBelowBed); a row
%         whose moments are then not finite - every trace blanked - is NOT
%         ice and ABSTAINS: theta, dlam, Cmag, phi, dP_hh and dP_hv are NaN
%         there, the orientation search never sees it, and the depth
%         kernels (coherence window, gradient) treat it as absent - zero
%         weight in every sum, and both observables are ratios in which
%         the window's finite-row count cancels - so the ice above keeps
%         its own value rather than inheriting NaN from below, with the
%         same one-sided-kernel edge behaviour the record's own top has
%         (dlam reads low within ~2 sigma of either). With no such row
%         every output is bit-identical to the unmasked chain.
%   z     [Nt x 1] depth of each row [m]
%   opts  fc (750e6), psi_step_deg (1), win_m (coherence window, 30),
%         grad_win_m (Gaussian sigma for the gradient, 25),
%         coh_min (0.4), deramped (true)
%
% Output fields: psi, dP_hh, dP_hv, phi, Cmag, psi_grad, dlam, theta,
%   theta_quality, plus the constants this call ran under - coh_min,
%   deramped, fc, eps_perp, deps, win_m. The last four are RECORDED rather
%   than merely defaulted because dlam is derived from them (eq. 9): any
%   downstream step that accepts dlam without re-fitting it, such as
%   ptt.ershadiInverse, has to propagate phase under the same constants,
%   and a mismatch there is a silent rescale with nothing to reveal it.

if nargin < 3, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
step = H_opt(opts, 'psi_step_deg', 1);
win_m = H_opt(opts, 'win_m', 30);
grad_win_m = H_opt(opts, 'grad_win_m', 25);
coh_min = H_opt(opts, 'coh_min', 0.4);
deramped = H_opt(opts, 'deramped', true);

z = z(:);
Nt = numel(z);
dz = median(abs(diff(z)));
if ~isfield(S, 'vh'), S.vh = S.hv; end   % reciprocity, as the paper assumes

C = ptt.constants();
eps_perp = 3.15;          % the paper's value, not eps_bar
deps = C.deps;            % 0.034, same as ptt.constants
c0 = C.c * 1e9;           % m/s
psi = (0:step:180-step) * pi/180;
Np = numel(psi);

% --- moments once; every observable below is quadratic in S (see
% ptt.quadpolMoments for why this is not done on the images)
M = ptt.quadpolMoments(S, [1 size(S.hh, 2)]);

Phh = zeros(Nt, Np); Pvv = zeros(Nt, Np); Phv = zeros(Nt, Np);
Chhvv = complex(zeros(Nt, Np));
for j = 1:Np
  % (E1) R S R' rather than R' S R: same algebra with psi -> -psi
  cc = cos(-psi(j)); ss = sin(-psi(j)); cs = cc*ss;
  w_hh = [cc*cc, ss*ss,  cs,  cs];
  w_vv = [ss*ss, cc*cc, -cs, -cs];
  w_hv = [-cs,     cs, cc*cc, -ss*ss];
  Phh(:, j) = real(H_quad(M, w_hh, w_hh));
  Pvv(:, j) = real(H_quad(M, w_vv, w_vv));
  Phv(:, j) = real(H_quad(M, w_hv, w_hv));
  Chhvv(:, j) = H_quad(M, w_hh, w_vv);
end

% --- rows with no finite moment (every trace blanked there) abstain from
% everything below; a row is either whole or absent
fin = all(isfinite(Phh) & isfinite(Pvv) & isfinite(Phv) & isfinite(Chhvv), 2);
Phh(~fin, :) = 0; Pvv(~fin, :) = 0; Phv(~fin, :) = 0; Chhvv(~fin, :) = 0;

% --- (E4) power anomalies, eq. (12): dB against the azimuthal mean of the
% AMPLITUDE, so the reference is mean|s| and not mean|s|^2
A_hh = sqrt(max(Phh, 0));
A_hv = sqrt(max(Phv, 0));
dP_hh = 20*log10(max(A_hh, realmin) ./ max(mean(A_hh, 2), realmin));
dP_hv = 20*log10(max(A_hv, realmin) ./ max(mean(A_hv, 2), realmin));
dP_hh(~fin, :) = NaN;
dP_hv(~fin, :) = NaN;

% --- coherence over a depth window, eq. (7), then (E2) the deramp conjugate.
% Absent rows are zero in numerator and denominator alike, so the ratio is
% the coherence over the finite rows in the window with no normaliser.
nw = max(3, 2*floor(win_m / max(dz, eps) / 2) + 1);
k = ones(nw, 1) / nw;
num = conv2(real(Chhvv), k, 'same') + 1i*conv2(imag(Chhvv), k, 'same');
den = sqrt(max(conv2(Phh, k, 'same'), 0) .* max(conv2(Pvv, k, 'same'), 0));
Cn = num ./ max(den, realmin);
if deramped
  Cn = conj(Cn);
end
Cn(~fin, :) = NaN;
Cmag = abs(Cn);
phi = angle(Cn);

% --- (E3) unwrap-free gradient: convolve the real and imaginary parts,
% then dphi/dz = Im(conj(C) dC/dz)/|C|^2. Absent rows enter both sums as
% zero, exactly as the rows beyond the record's ends always have; a
% normaliser over the finite rows would be a real positive factor common
% to Cs and dC and drops out of the ratio, so none is applied.
ns = max(3, 2*ceil(3 * grad_win_m / max(dz, eps)) + 1);
x = (-(ns-1)/2:(ns-1)/2).' * dz;
sig = max(grad_win_m, dz);
g = exp(-0.5*(x/sig).^2); g = g / sum(g);
dg = -(x./sig.^2) .* g; dg = dg - mean(dg);      % zero-sum derivative kernel
Cz = Cn; Cz(~fin, :) = 0;
Cs = conv2(real(Cz), g, 'same') + 1i*conv2(imag(Cz), g, 'same');
dC = conv2(real(Cz), dg, 'same') + 1i*conv2(imag(Cz), dg, 'same');
dphi = imag(conj(Cs) .* dC) ./ max(abs(Cs).^2, realmin);
dphi(~fin, :) = NaN;

% --- eq. (9)/(10) with the square-root form
psi_grad = (c0 * sqrt(eps_perp) / (2*pi*fc*deps)) * dphi;

% --- (E5) quality gate.
% The gate is DESTRUCTIVE and it is applied to the direct chain's own
% products only. Cn is kept ungated in the output as C_raw, because the
% strength stage (ptt.ershadiStrength) fits the complex coherence WEIGHTED
% by |C| rather than thresholded on it: a hard cut at 0.4 removes rows that
% still carry usable phase - at Ridge A, whose median |C| is 0.47, it voids
% 61% of deep intervals and leaves the strength profile with nothing to be
% inferred from. Weighting degrades where the data degrades; a gate falls
% off a cliff.
C_raw = Cn;
bad = Cmag < coh_min;
psi_grad(bad) = NaN;
phi(bad) = NaN;

% --- orientation: minimum of dP_hv, which has period 90 deg, then the
% polarity of phi at that azimuth to choose between the two axes
theta = nan(Nt, 1);
[~, imin] = min(dP_hv(fin, :), [], 2);
% NEGATED. The sweep is built with the paper's R S R' sense, which is the
% inverse rotation to ptt.rotatePolarization's R' S R, so the sweep index
% that aligns with the fabric sits at MINUS the physical azimuth. Reporting
% the index angle directly returned 74.7 deg for a synthetic truth of 15
% and 20.3 for a truth of 70, while every apparent pass was at an angle
% that is its own negative modulo 90 (0 and 45) - which is exactly how a
% sign convention hides until it is tested off the symmetry points. The
% printed equation cannot settle this on its own, because it turns on how
% R is defined; the synthetic does.
theta(fin) = mod(-psi(imin).', pi);
theta_quality = max(dP_hv, [], 2) - min(dP_hv, [], 2);   % dB modulation
dl_at = nan(Nt, 1);
dl_at(fin) = psi_grad(sub2ind([Nt Np], find(fin), imin));
% Eq. (10) evaluates Psi at the principal axes; the axis carrying lam_max
% is the one where Psi is positive, so a negative value means the minimum
% found was the other axis and theta moves by 90 deg.
flip = dl_at < 0;
theta(flip) = mod(theta(flip) + pi/2, pi);
dlam = abs(dl_at);

out = struct('psi', psi, 'dP_hh', dP_hh, 'dP_hv', dP_hv, 'phi', phi, ...
  'Cmag', Cmag, 'C_raw', C_raw, 'psi_grad', psi_grad, 'dlam', dlam, 'theta', theta(:), ...
  'theta_quality', theta_quality, 'coh_min', coh_min, ...
  'deramped', deramped, 'fc', fc, 'eps_perp', eps_perp, 'deps', deps, ...
  'win_m', win_m);

end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function v = H_quad(M, w1, w2)
v = complex(zeros(size(M,1), 1));
for k = 1:4
  if w1(k) == 0, continue; end
  for l = 1:4
    if w2(l) == 0, continue; end
    v = v + (w1(k) * w2(l)) * M(:, k, l);
  end
end
end
