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
%   S     struct of complex [Nt x Nx] channels hh, vv, hv (vh optional)
%   z     [Nt x 1] depth of each row [m]
%   opts  fc (750e6), psi_step_deg (1), win_m (coherence window, 30),
%         grad_win_m (Gaussian sigma for the gradient, 25),
%         coh_min (0.4), deramped (true)
%
% Output fields: psi, dP_hh, dP_hv, phi, Cmag, psi_grad, dlam, theta,
%   theta_quality

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

% --- (E4) power anomalies, eq. (12): dB against the azimuthal mean of the
% AMPLITUDE, so the reference is mean|s| and not mean|s|^2
A_hh = sqrt(max(Phh, 0));
A_hv = sqrt(max(Phv, 0));
dP_hh = 20*log10(max(A_hh, realmin) ./ max(mean(A_hh, 2), realmin));
dP_hv = 20*log10(max(A_hv, realmin) ./ max(mean(A_hv, 2), realmin));

% --- coherence over a depth window, eq. (7), then (E2) the deramp conjugate
nw = max(3, 2*floor(win_m / max(dz, eps) / 2) + 1);
k = ones(nw, 1) / nw;
num = conv2(real(Chhvv), k, 'same') + 1i*conv2(imag(Chhvv), k, 'same');
den = sqrt(max(conv2(Phh, k, 'same'), 0) .* max(conv2(Pvv, k, 'same'), 0));
Cn = num ./ max(den, realmin);
if deramped
  Cn = conj(Cn);
end
Cmag = abs(Cn);
phi = angle(Cn);

% --- (E3) unwrap-free gradient: convolve the real and imaginary parts,
% then dphi/dz = Im(conj(C) dC/dz)/|C|^2
ns = max(3, 2*ceil(3 * grad_win_m / max(dz, eps)) + 1);
x = (-(ns-1)/2:(ns-1)/2).' * dz;
sig = max(grad_win_m, dz);
g = exp(-0.5*(x/sig).^2); g = g / sum(g);
dg = -(x./sig.^2) .* g; dg = dg - mean(dg);      % zero-sum derivative kernel
Cs = conv2(real(Cn), g, 'same') + 1i*conv2(imag(Cn), g, 'same');
dC = conv2(real(Cn), dg, 'same') + 1i*conv2(imag(Cn), dg, 'same');
dphi = imag(conj(Cs) .* dC) ./ max(abs(Cs).^2, realmin);

% --- eq. (9)/(10) with the square-root form
psi_grad = (c0 * sqrt(eps_perp) / (2*pi*fc*deps)) * dphi;

% --- (E5) quality gate
bad = Cmag < coh_min;
psi_grad(bad) = NaN;
phi(bad) = NaN;

% --- orientation: minimum of dP_hv, which has period 90 deg, then the
% polarity of phi at that azimuth to choose between the two axes
[~, imin] = min(dP_hv, [], 2);
% NEGATED. The sweep is built with the paper's R S R' sense, which is the
% inverse rotation to ptt.rotatePolarization's R' S R, so the sweep index
% that aligns with the fabric sits at MINUS the physical azimuth. Reporting
% the index angle directly returned 74.7 deg for a synthetic truth of 15
% and 20.3 for a truth of 70, while every apparent pass was at an angle
% that is its own negative modulo 90 (0 and 45) - which is exactly how a
% sign convention hides until it is tested off the symmetry points. The
% printed equation cannot settle this on its own, because it turns on how
% R is defined; the synthetic does.
theta = mod(-psi(imin).', pi);
theta_quality = max(dP_hv, [], 2) - min(dP_hv, [], 2);   % dB modulation
lin = sub2ind([Nt Np], (1:Nt).', imin);
dl_at = psi_grad(lin);
% Eq. (10) evaluates Psi at the principal axes; the axis carrying lam_max
% is the one where Psi is positive, so a negative value means the minimum
% found was the other axis and theta moves by 90 deg.
flip = dl_at < 0;
theta(flip) = mod(theta(flip) + pi/2, pi);
dlam = abs(dl_at);

out = struct('psi', psi, 'dP_hh', dP_hh, 'dP_hv', dP_hv, 'phi', phi, ...
  'Cmag', Cmag, 'psi_grad', psi_grad, 'dlam', dlam, 'theta', theta(:), ...
  'theta_quality', theta_quality, 'coh_min', coh_min, ...
  'deramped', deramped);

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
