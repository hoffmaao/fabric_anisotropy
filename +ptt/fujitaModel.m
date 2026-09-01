function out = fujitaModel(layers, z, psi, opts)
%FUJITAMODEL Layered radar forward model of Fujita et al. (2006).
%
% out = ptt.fujitaModel(layers, z, psi, opts)
%
% The matrix model of Fujita, Maeno and Matsuoka (2006, J. Glaciol. 52,
% 407-424), in the form Ershadi et al. (2022, eq. 5) fit it:
%
%   S_N = D^2 [prod R T R']_down  R(th_N) Gamma R'(th_N)  [prod R T R']_up
%
% evaluated for every scattering depth z and every antenna azimuth psi.
% Each layer i contributes a one-way transmission T_i = diag(e^{i k_x dz},
% e^{i k_y dz}) in its own principal frame at th_i, and the boundary at the
% bottom of the stack scatters with Gamma = diag(Gamma_x, Gamma_y), also in
% the local frame. Everything common to both polarizations - the depth
% factor D^2, absolute attenuation, the mean propagation phase - drops out
% of every observable this returns (power anomalies are normalized over
% azimuth, the coherence phase is relative), so only the RELATIVE phase
% delta_i = gpd*dlam_i*dz and the RATIO r_i = Gamma_y/Gamma_x enter.
%
% CONVENTIONS, chosen to match the fields ptt.ershadiFabric extracts from
% data so a model evaluated here is directly comparable:
%   - psi is the SWEEP azimuth of the synthesized antennas, applied in the
%     paper's R S R' sense (ershadiFabric E1): a fabric at physical
%     azimuth theta produces features at sweep index -theta. This file
%     takes the same sense, so model and data line up index for index.
%   - The coherence is eq. (7) over the same depth window win_m the data
%     path uses, and it is NOT conjugated: the paper uses eq. (7) for
%     models and its conjugate for radar data (ershadiFabric E2 applies
%     that conjugate to the data), so both end in the same sense.
%   - r_dB = 20*log10(Gamma_y/Gamma_x). Amplitude-coefficient convention,
%     fixed by their eq. (13): at the anti-phase depth the co-pol nodes of
%     a uniform layer sit at tan^2(AD/2) = Gamma_x/Gamma_y, and their
%     Fig. 4 curve reproduces under 20log10, not 10log10.
%
% Inputs
%   layers  struct array, one per layer, TOP FIRST:
%             .top_m   top depth of the layer (first must be 0)
%             .dlam    horizontal eigenvalue difference in the layer
%             .theta   physical axis azimuth, radians, mod pi
%             .r_db    20log10(Gamma_y/Gamma_x) at the layer's bottom
%                      boundary (and at every depth within it)
%   z       depth axis (m), column
%   psi     sweep azimuths (radians, row or column)
%   opts    .fc (750e6), .eps_perp (3.15, the paper's value),
%           .deps (0.034), .win_m (30, eq.-7 window, match the data path)
%
% Output fields (Nz x Npsi unless noted)
%   s_hh, s_vv, s_hv   complex scattering amplitudes (relative units)
%   dP_hh, dP_hv       eq. (12) power anomalies, dB
%   C                  eq. (7) coherence over win_m
%   phi                angle(C)
%   Cmag               abs(C)
%
% See also ptt.ershadiFabric, ptt.ershadiInverse.

if nargin < 4, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
eps_perp = H_opt(opts, 'eps_perp', 3.15);
deps = H_opt(opts, 'deps', 0.034);
win_m = H_opt(opts, 'win_m', 30);

C0 = ptt.constants();
c0 = C0.c * 1e9;
% ONE-WAY relative phase per metre per unit dlam; the down+up sandwich
% doubles it, recovering the two-way gpd the rest of the toolbox uses.
gpd1 = pi * fc * deps / (sqrt(eps_perp) * c0);

z = z(:);
Nz = numel(z);
psi = psi(:).';
Np = numel(psi);
tops = [layers.top_m];
if tops(1) ~= 0
  error('ptt:fujitaModel:layers', 'first layer must start at 0 m');
end
if any(diff(tops) <= 0)
  error('ptt:fujitaModel:layers', 'layer tops must strictly increase');
end
NL = numel(layers);

% --- cumulative one-way Jones matrix of the stack above each depth, in the
% ANTENNA frame, bottom boundary of the stack scattering with the local
% Gamma. Layers are uniform, so within layer L at depth z the product is
%   P(z) = A_L(z - top_L) * A_{L-1}(dz_{L-1}) * ... * A_1(dz_1)
% with A_i(d) = R(th_i) diag(e^{+i gpd1 dlam_i d/2}, e^{-i ...}) R(th_i)'.
% The scattered return is S = P.' * G * P (each A_i is symmetric, so the
% up-path product is the transpose of the down one - eq. (5) collapses to
% this for reciprocal media).
s_hh = complex(zeros(Nz, Np));
s_vv = complex(zeros(Nz, Np));
s_hv = complex(zeros(Nz, Np));

% per-depth layer index and boundary products
Pacc = repmat(eye(2), 1, 1);            % product of the FULL layers above
li = 1;
for iz = 1:Nz
  while li < NL && z(iz) >= tops(li+1)
    % close layer li: multiply its full-thickness A into the accumulator
    d = tops(li+1) - tops(li);
    Pacc = H_A(layers(li), d, gpd1) * Pacc;
    li = li + 1;
  end
  d = z(iz) - tops(li);
  P = H_A(layers(li), d, gpd1) * Pacc;
  r_amp = 10^(layers(li).r_db / 20);
  G = H_R(layers(li).theta) * diag([1, r_amp]) * H_R(layers(li).theta).';
  M = P.' * G * P;                       % physical two-way Jones at z
  % sweep in the paper's R S R' sense (ershadiFabric E1): psi -> -psi
  cc = cos(-psi); ss = sin(-psi);
  m11 = M(1,1); m12 = M(1,2); m21 = M(2,1); m22 = M(2,2);
  s_hh(iz,:) = cc.*cc*m11 + ss.*ss*m22 + cc.*ss*(m12 + m21);
  s_vv(iz,:) = ss.*ss*m11 + cc.*cc*m22 - cc.*ss*(m12 + m21);
  s_hv(iz,:) = -cc.*ss*m11 + cc.*ss*m22 + cc.*cc*m12 - ss.*ss*m21;
end

% --- eq. (12) power anomalies against the azimuthal mean AMPLITUDE
A_hh = abs(s_hh); A_hv = abs(s_hv);
out.dP_hh = 20*log10(max(A_hh, realmin) ./ max(mean(A_hh, 2), realmin));
out.dP_hv = 20*log10(max(A_hv, realmin) ./ max(mean(A_hv, 2), realmin));

% --- eq. (7) coherence over the same depth window as the data path; the
% model is deterministic so the window only mimics the data's smoothing
dz = median(abs(diff(z)));
nw = max(3, 2*floor(win_m / max(dz, eps) / 2) + 1);
k = ones(nw, 1) / nw;
num = conv2(real(s_hh .* conj(s_vv)), k, 'same') + ...
  1i*conv2(imag(s_hh .* conj(s_vv)), k, 'same');
den = sqrt(max(conv2(abs(s_hh).^2, k, 'same'), 0) .* ...
  max(conv2(abs(s_vv).^2, k, 'same'), 0));
Cn = num ./ max(den, realmin);
out.C = Cn;
out.phi = angle(Cn);
out.Cmag = abs(Cn);
out.s_hh = s_hh; out.s_vv = s_vv; out.s_hv = s_hv;
out.psi = psi;

end

function A = H_A(L, d, gpd1)
%H_A One-way Jones of thickness d of layer L, antenna frame. Symmetric.
del = gpd1 * L.dlam * d;                 % one-way x-vs-y relative phase
R = H_R(L.theta);
A = R * diag([exp(+0.5i*del), exp(-0.5i*del)]) * R.';
end

function R = H_R(t)
R = [cos(t), -sin(t); sin(t), cos(t)];
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
