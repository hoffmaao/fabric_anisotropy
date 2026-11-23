function [S, tau, alpha] = synth_quadpol_column(z, dlam, theta, dt, p)
%SYNTH_QUADPOL_COLUMN Four antenna-frame channels of a birefringent column with speckle.
%
% [S, tau, alpha] = synth_quadpol_column(z, dlam, theta, dt, p)
%
% A slow-axis signal a and a fast-axis signal b of the SAME reflectivity,
% a delayed by tau(z) = k * int dlam dz with the carrier phase
% exp(-2 pi i fc tau) that mixing by the LO leaves on a delayed echo (the
% ENGINEERING sign, as OPR stores data), then rotated into an antenna
% frame whose slow axis sits at theta [rad] (ptt.rotatePolarization's
% sense: rotating S by theta returns hh = a, vv = b):
%   hh = a c^2 + b s^2,  vv = a s^2 + b c^2,  hv = vh = (a - b) c s.
% Then channel gains (receive gr on vv and hv, transmit gt on vv and vh, hh
% the reference), and white noise at p.snr_db at the top falling 15 dB to
% the bottom.
%
% DEPOLARISED BACKSCATTER (p.depol > 0). Below the firn real ice does not
% return one polarised echo per depth: the scatterer is a random symmetric
% matrix in the fabric's own frame,
%   [c + d3, d1; d1, c - d3],
% c the coherent part above and d3, d1 a depolarised part, independent of
% it and of each other, each of relative power p.depol. Every element goes
% through the same birefringent column, so in the fabric frame
%   a = A[c] + A[d3]   (slow-slow: delayed by tau)
%   b = B[c] - B[d3]   (fast-fast: undelayed)
%   x = H[d1]          (cross: one leg each, delayed by tau / 2)
% and the antenna frame is the rotation of [a x; x b]. The eigen-pair
% coherence falls to (1 - depol) / (1 + depol) and the cross-polarised
% power in the fabric frame rises to depol / (1 + depol) of the
% co-polarised, with NO bias on the pair's phase: the depolarised echo
% carries the column's delay like the coherent one. (Adding independent
% fields to the antenna-frame channels instead, as this helper once did,
% leaves a zero-phase term in any rotated frame and a phase error that
% follows the fringes, which crossing passes over real ice do not show.)
% depol 0.4 gives the coherence (0.43) and the in-frame cross-polarised
% power (-5 dB) of Ridge A below the firn.
%
% p fields (defaults): fc 750e6, f_off -13e6 (band offset from the LO as
% stored), B 160e6, Nx 125, depol 0, snr_db 15, tau_x 3 (speckle
% correlation, traces), gr 1, gt 1, refl_phase_deg 0, refl_len_m 30.
% z must be uniform in fast time: z = (0:Nt-1) * dz.
%
% The band offset matters for the nodes: a stored baseband frequency f is
% the radio frequency fc + f, so the cross-polarised power 45 deg off the
% axes turns once per fringe of the band's CENTRE, fc + f_off, while the
% carrier phase read at the envelope peak counts fringes of fc.
%
% refl_phase_deg > 0 gives the slow-axis reflectivity a phase against the
% fast-axis one, the same for every trace and varying smoothly with depth
% (standard deviation refl_phase_deg, its autocorrelation falling to 1/e
% at refl_len_m): the
% layers' own contribution to the phase between the polarisations. It does
% not accumulate - it is there at each depth and gone below it - which is
% what separates it from a delay. [S, tau, alpha] returns it in radians.
fc = H_opt(p, 'fc', 750e6); f_off = H_opt(p, 'f_off', -13e6); B = H_opt(p, 'B', 160e6);
Nx = H_opt(p, 'Nx', 125); depol = H_opt(p, 'depol', 0); snr_db = H_opt(p, 'snr_db', 15);
tau_x = H_opt(p, 'tau_x', 3); gr = H_opt(p, 'gr', 1); gt = H_opt(p, 'gt', 1);
refl_deg = H_opt(p, 'refl_phase_deg', 0); refl_len = H_opt(p, 'refl_len_m', 30);
C = ptt.constants();
k = C.deps / (sqrt(C.eps_bar) * C.c);
z = z(:); dlam = dlam(:);
tau = k * cumtrapz(z, dlam);
Nt = numel(z);
Nf = 2 ^ nextpow2(2 * Nt);
f = ifftshift(((0:Nf-1) - floor(Nf / 2)) / (Nf * dt)).';     % GHz
fo = f - f_off / 1e9;
pulse = double(abs(fo) <= B / 2e9) .* (0.5 + 0.5 * cos(pi * fo / (B / 2e9)));
amp = 10 .^ (-(z / z(end)) * 15 / 20);
r = H_speckle(Nt, Nx, tau_x) .* amp;
% the reflectors' own phase between the two polarisations: white noise
% smoothed by a Gaussian of half the correlation length (its
% autocorrelation is then exp(-(lag / refl_len)^2)) and scaled to the
% standard deviation
alpha = zeros(Nt, 1);
if refl_deg > 0
  n_c = max(1, round(refl_len / 2 / median(diff(z))));
  kk = exp(-(-3*n_c:3*n_c).' .^ 2 / (2 * n_c ^ 2));
  alpha = conv(randn(Nt + 6 * n_c, 1), kk, 'valid');
  alpha = deg2rad(refl_deg) * alpha / std(alpha);
end
% the coherent echo, and with it the reflectors' phase on the slow axis
[a, b] = H_through(r, tau, alpha, f, fc, pulse, Nf);
x = zeros(Nt, Nx);
if depol > 0
  % the depolarised part of the scatterer, through the same column
  [ad, bd] = H_through(H_speckle(Nt, Nx, tau_x) .* amp * sqrt(depol), tau, [], f, fc, pulse, Nf);
  a = a + ad; b = b - bd;
  x = H_through(H_speckle(Nt, Nx, tau_x) .* amp * sqrt(depol), tau / 2, [], f, fc, pulse, Nf);
end
% the fabric frame [a x; x b] seen from antennas whose slow axis sits at theta
c = cos(theta); s = sin(theta);
S.hh = a * c ^ 2 + b * s ^ 2 - 2 * c * s * x;
S.vv = a * s ^ 2 + b * c ^ 2 + 2 * c * s * x;
S.hv = (a - b) * c * s + (c ^ 2 - s ^ 2) * x;
S.vh = S.hv;
S.vv = S.vv * gr * gt; S.hv = S.hv * gr; S.vh = S.vh * gt;
nz = sqrt(mean(abs(b(1:min(200, Nt), :)) .^ 2, 'all')) * 10 ^ (-snr_db / 20);
for ch = {'hh', 'vv', 'hv', 'vh'}
  S.(ch{1}) = S.(ch{1}) + nz * (randn(Nt, Nx) + 1i * randn(Nt, Nx)) / sqrt(2);
end
end

function [a, b] = H_through(r, tau, alpha, f, fc, pulse, Nf)
% One reflectivity field r [Nt x Nx] as two band-limited echoes: b
% undelayed, a delayed by tau(z), built in chunks over which tau is
% constant, each chunk shifted in the frequency domain and given the phase
% exp(-2 pi i fc tau) that mixing by the LO leaves on a delayed echo. alpha
% (radians, [Nt x 1], or empty) is a phase on a's reflectivity alone.
[Nt, Nx] = size(r);
chunk = 40;
a = zeros(Nt, Nx); b = zeros(Nt, Nx);
for s0 = 1:chunk:Nt
  i = s0:min(Nt, s0 + chunk - 1);
  ti = mean(tau(i));
  ri = zeros(Nt, Nx); ri(i, :) = r(i, :);
  R = fft(ri, Nf, 1) .* pulse;
  if nargout > 1, b = b + H_crop(ifft(R, [], 1), Nt); end
  if ~isempty(alpha) && any(alpha(i))
    ri(i, :) = r(i, :) .* exp(-1i * alpha(i));
    R = fft(ri, Nf, 1) .* pulse;
  end
  a = a + H_crop(ifft(R .* exp(-2i * pi * f * ti), [], 1), Nt) * exp(-2i * pi * fc * 1e-9 * ti);
end
end

function x = H_crop(x, Nt)
x = x(1:Nt, :);
end

function s = H_speckle(Nt, Nx, tau_x)
s = (randn(Nt, Nx) + 1i * randn(Nt, Nx)) / sqrt(2);
if tau_x > 1
  kk = exp(-(-3*tau_x:3*tau_x) .^ 2 / (2 * (tau_x / 2) ^ 2)); kk = kk / norm(kk);
  s = conv2(1, kk, s, 'same');
end
end

function v = H_opt(o, name, default)
if isfield(o, name) && ~isempty(o.(name)), v = o.(name); else, v = default; end
end
