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
%             .gx_db   OPTIONAL amplitude of Gamma_x, dB re the unit
%                      reference (default 0): the layer's absolute
%                      scattering strength. dP_hh and dP_hv are invariant
%                      to it exactly, being normalized over azimuth WITHIN
%                      a row. The eq.-(7) coherence is invariant within a
%                      row for the same reason (it is a ratio), but its
%                      win_m depth window MIXES rows, so a layered gx_db
%                      profile does reweight C over the ~win_m straddling
%                      each layer boundary - it is not invariant there.
%                      gx_db exists for the power PROFILE below, which is
%                      where relative reflection strengths live.
%   z       depth axis (m), column
%   psi     sweep azimuths (radians, row or column)
%   opts    .fc (750e6), .eps_perp (3.15, the paper's value),
%           .deps (0.034), .win_m (30, eq.-7 window, match the data path)
%           .bed OPTIONAL struct('z_m', depth, 'gx_db', amp, 'r_db', 0):
%                a single strong interface terminating the domain - the
%                ice-bed reflection, tens of dB above the internals
%                (Fujita's Table-2 internals sit at Gamma_x = 1e-12;
%                a bed amplitude coefficient of 0.1-0.3 is typical). The
%                row CONTAINING bed z_m scatters with the bed's Gamma in
%                the local layer frame; every row BELOW it returns NaN in
%                all fields, because the domain ends at the bed - sub-bed
%                returns are the thing the movies grey out, not a
%                prediction this model should make. The eq.-(7) coherence
%                is windowed over the rows ABOVE the bed only, so the NaN
%                band is exactly the sub-bed rows and does not bleed
%                (nw-1)/2 rows upward; the bed's own row and the ~win_m/2
%                above it are still window EDGE values, damped by the
%                zero padding at the truncation exactly as the bottom of
%                a bedless depth axis is.
%
% Output fields (Nz x Npsi unless noted)
%   s_hh, s_vv, s_hv   complex scattering amplitudes (relative units)
%   dP_hh, dP_hv       eq. (12) power anomalies, dB
%   C                  eq. (7) coherence over win_m
%   phi                angle(C)
%   Cmag               abs(C)
%   P_hh_db            azimuthal-mean RELATIVE power profile (Nz x 1),
%                      10log10(mean_psi |s_hh|^2): internals sit at their
%                      gx_db and the bed spikes above them - the relative
%                      strength of the reflections, which no normalized
%                      observable can show
%   gpd1               the one-way phase rate this evaluation used (scalar)
%
% See also ptt.ershadiFabric, ptt.ershadiInverse.

if nargin < 4, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
eps_perp = H_opt(opts, 'eps_perp', 3.15);
deps = H_opt(opts, 'deps', 0.034);
win_m = H_opt(opts, 'win_m', 30);

% ONE-WAY relative phase per metre per unit dlam; the down+up sandwich
% doubles it, recovering the two-way gpd the rest of the toolbox uses.
% ptt.ershadiInverse locates its eq.-(13) anti-phase depths with the same
% call, so the two cannot drift apart.
gpd1 = ptt.birefringentPhaseRate(fc, eps_perp, deps);

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

% optional bed: a single strong interface; the row containing it takes the
% bed's Gamma, rows below it have no return at all
bed = H_opt(opts, 'bed', []);
bed_row = 0;
if ~isempty(bed)
  if ~isfield(bed, 'z_m') || ~isfield(bed, 'gx_db')
    error('ptt:fujitaModel:bed', 'opts.bed needs z_m and gx_db');
  end
  if ~isfield(bed, 'r_db') || isempty(bed.r_db), bed.r_db = 0; end
  bed_row = find(z >= bed.z_m, 1);
  if isempty(bed_row), bed_row = 0; end   % bed below the axis: no effect
end

% per-depth layer index and boundary products
Pacc = repmat(eye(2), 1, 1);            % product of the FULL layers above
li = 1;
for iz = 1:Nz
  if bed_row > 0 && iz > bed_row
    s_hh(iz,:) = NaN; s_vv(iz,:) = NaN; s_hv(iz,:) = NaN;
    continue                             % the domain ends at the bed
  end
  while li < NL && z(iz) >= tops(li+1)
    % close layer li: multiply its full-thickness A into the accumulator
    d = tops(li+1) - tops(li);
    Pacc = H_A(layers(li), d, gpd1) * Pacc;
    li = li + 1;
  end
  d = z(iz) - tops(li);
  P = H_A(layers(li), d, gpd1) * Pacc;
  if bed_row > 0 && iz == bed_row
    gx = 10^(bed.gx_db / 20);
    r_amp = 10^(bed.r_db / 20);
  else
    gx = 10^(H_gx(layers(li)) / 20);
    r_amp = 10^(layers(li).r_db / 20);
  end
  G = H_R(layers(li).theta) * (gx * diag([1, r_amp])) * H_R(layers(li).theta).';
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
% model is deterministic so the window only mimics the data's smoothing.
% THE WINDOW IS FORMED ABOVE THE BED, NOT ACROSS IT. conv2 propagates a NaN
% over the whole kernel reach, so running it on a column whose sub-bed rows
% are NaN would blank C, phi and Cmag for the (nw-1)/2 rows ABOVE the bed
% too - a 15 m band at the default win_m = 30 that a caller would read as
% decoherence rather than as the model's own boundary handling. Truncating
% at bed_row instead leaves every row more than (nw-1)/2 above the bed
% bit-identical (the kernel never reaches the cut) and gives the rows next
% to the bed the same zero-padded edge treatment the bottom of a bedless
% depth axis already gets.
dz = median(abs(diff(z)));
nw = max(3, 2*floor(win_m / max(dz, eps) / 2) + 1);
k = ones(nw, 1) / nw;
if bed_row > 0, iw = 1:bed_row; else, iw = 1:Nz; end
sh = s_hh(iw, :); sv = s_vv(iw, :);
num = conv2(real(sh .* conj(sv)), k, 'same') + ...
  1i*conv2(imag(sh .* conj(sv)), k, 'same');
den = sqrt(max(conv2(abs(sh).^2, k, 'same'), 0) .* ...
  max(conv2(abs(sv).^2, k, 'same'), 0));
Cn = nan(Nz, Np);
Cn(iw, :) = num ./ max(den, realmin);
out.C = Cn;
out.phi = angle(Cn);
out.Cmag = abs(Cn);
out.s_hh = s_hh; out.s_vv = s_vv; out.s_hv = s_hv;
out.psi = psi;
% relative reflection strength vs depth: internals at their gx_db, the bed
% spiking above them. 10log10 of mean power = 20log10 of amplitude scale,
% so a +30 dB gx_db bed reads ~+30 dB here.
out.P_hh_db = 10*log10(max(mean(abs(s_hh).^2, 2), realmin));
if bed_row > 0
  out.P_hh_db(bed_row+1:end) = NaN;
end
out.gpd1 = gpd1;

end

function g = H_gx(L)
if isfield(L, 'gx_db') && ~isempty(L.gx_db), g = L.gx_db; else, g = 0; end
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
