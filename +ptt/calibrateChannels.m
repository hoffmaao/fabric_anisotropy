function [Mc, g, info] = calibrateChannels(M, z, opts)
%CALIBRATECHANNELS Complex per-channel gains of a quad-pol system, from its own moments.
%
% [Mc, g, info] = ptt.calibrateChannels(M, z, opts)
%
% Removes the transmit and receive gains of the V channel relative to the H
% channel - amplitude AND phase - from an antenna-frame moment matrix, so
% that every azimuth synthesised from it is a combination of channels on
% one scale. With a = t_v/t_h (transmit) and b = r_v/r_h (receive), the
% measured channels are
%
%     HH = S_hh,   VV = a b S_vv,   HV = a S_hv,   VH = b S_vh
%
% (HV: transmit V, receive H), so the measured moments are
% M_kl = g_k conj(g_l) M_true_kl with g = [1, ab, a, b], and the
% calibration is a division of each moment by g_k conj(g_l). Four real
% numbers are unknown - |a|, |b|, arg a, arg b - and the ice supplies
% three constraints and the firn a fourth:
%
%   RECIPROCITY. A reciprocal medium and antennas give S_hv = S_vh, so
%   <HV VH*> = a conj(b) <|S_hv|^2>: its phase is arg a - arg b and the
%   power ratio |HV|^2/|VH|^2 is |a|^2/|b|^2, at every depth. Measured on
%   every season this system has flown, arg a - arg b is +110 to +112
%   deg (Ridge A, Thwaites/WAIS/McMurdo, EastGRIP), with the cross-pol
%   pair coherent at 0.94-0.99, so this is a design constant of the
%   instrument, not a per-frame accident.
%
%   THE FIRN, AMPLITUDE. In the firn window (default 60-150 m) the fabric
%   is weak and the reflections nearly isotropic, so <|VV|^2>/<|HH|^2>
%   should be 1 and what is measured is |ab|^2: -3.3 to -4.0 dB on every
%   Ridge A frame whatever its heading. With |b/a| from reciprocity that
%   fixes |a| and |b| separately.
%
%   THE FIRN, PHASE. The co-pol phase arg<HH VV*> at depth z is
%   -arg(ab) plus the birefringent phase, which is zero at the surface
%   and grows with depth at a rate set by the fabric and the heading. Its
%   value at z = 0 is therefore -arg(ab) alone. It is read by fitting the
%   single-axis form of that phase through the firn window -
%   phi0 + atan2(mu sin d, A + B cos d), d = rate * z, mu = cos 2(psi - t)
%   the axis-to-antenna term with A = (1-mu^2)/2, B = (1+mu^2)/2 (the
%   same closed form ptt.quadpolFabricLS fits) - on a grid in (mu, rate)
%   with phi0 in closed form, because a straight line through the window
%   is only right while the phase there is small: a column carrying 0.06
%   from the surface has already turned 1-3 rad by 60-150 m, and the
%   line's intercept was 8 deg off on it. On Ridge A the phase at
%   30-60 m reads +55..+64 deg on the N-S lines and +26..+35 on the rows,
%   growing in opposite senses with depth as a birefringent term must,
%   about a common intercept; with three constraints in hand this is the
%   fourth and the only one that carries a fabric term, so it is the
%   least certain (see info.phase_fit) and the one the crossing-pair
%   test is there to check.
%
% THE SIGN OF a AND b. The four constraints fix ab and a/b, and so a and
% b only up to a COMMON sign: (-a, -b) has the same ab and a/b. The two
% branches are not equivalent. Negating both negates HV and VH against the
% co-pol channels, which mirrors every synthesised azimuth (psi -> -psi)
% and with it every antenna-frame axis, and no moment can tell them apart
% because a mirrored column is also a column. The branch taken is the one
% in which the V chains are, on average, in phase with the H chain -
% Re(a/|a| + b/|b|) >= 0 - which keeps the handedness the uncalibrated
% synthesis has always used, the one the Ridge A heading test validated
% (one geographic axis across five heading families). For this system
% (arg a - arg b = +111, arg ab = -45 deg) that is arg a = +33, arg b = -78
% rather than -147 and +102. A wrong branch would show as a mirrored axis
% between heading families, which the heading test would catch.
%
% WHY THIS MATTERS. Uncalibrated, the -3.6 dB cross-pol pedestal that
% antenna-locks every power-based axis estimate on this system is mostly
% these gains: on the same Ridge A moments the synthesised cross-pol null
% sits on the antennas raw (psi 0-6 deg on both heading families) and
% within 3-9 deg of the fabric axis calibrated, the LS pedestal a3 falls
% from ~0.27 to ~0, and the heading-family gap in dlam (N-S minus rows
% +0.0059, +4.4 deg in axis) closes to -0.0001 / +1.9 deg at a co-pol
% phase near 0 (issue #23). The 2023 Antarctic season and EastGRIP carry
% 10-30 dB of range-dependent amplitude imbalance on top (issue #24),
% which this does not model: it fits ONE complex gain per channel over
% the record. ptt.equaliseChannels handles range-dependent AMPLITUDE and
% is the right tool there; this handles the phases equaliseChannels
% cannot.
%
% Inputs
%   M     [Nt x 4 x 4] antenna-frame moments over (hh, vv, hv, vh), from
%         ptt.quadpolMoments, pooled over as many traces as the frame
%         has (the estimate is a survey constant; noise averages down)
%   z     [Nt x 1] depth [m]
%   opts  firn_m ([60 150]) the window for the amplitude anchor and the
%           phase intercept; recip_m ([150 1400]) the depth range the
%           reciprocity terms are summed over; fc (750e6) for the
%           birefringent phase rate the firn phase fit uses;
%         phase_ab ([] = fit the firn intercept; a scalar in RADIANS
%           pins arg(ab) to a value settled elsewhere, e.g. on the
%           crossing pairs, and the intercept is still reported)
%
% Outputs
%   Mc    the calibrated moments, M_kl / (g_k conj(g_l))
%   g     [1 x 4] complex gains [1, ab, a, b] that were divided out, on
%         the branch THE SIGN OF a AND b names
%   info  a_db, b_db (20log10 of the moduli), arg_a_deg, arg_b_deg,
%         xpol_phase_deg (arg a - arg b, from reciprocity), xpol_coh (the
%         HV-VH coherence it was read from), firn_db (10log10 of the
%         firn VV/HH power ratio), phase_fit (struct: z, phase_deg, the
%         fit's intercept_deg, its dlam and mu, its rms_deg, and the
%         straight line's intercept_deg for comparison), and
%         phase_ab_deg as applied, wrapped to (-180, 180]
%
% See also ptt.equaliseChannels, ptt.quadpolMoments, ptt.quadpolFabricLS.

if nargin < 3, opts = struct(); end
firn = H_opt(opts, 'firn_m', [60 150]);
recip = H_opt(opts, 'recip_m', [150 1400]);
phase_ab = H_opt(opts, 'phase_ab', []);
fc = H_opt(opts, 'fc', 750e6);
Cc = ptt.constants();
gpd = 2 * pi * fc * Cc.deps / (sqrt(Cc.eps_bar) * Cc.c * 1e9);   % rad/m per unit dlam, as the LS
z = z(:);
if ndims(M) ~= 3 || size(M, 2) ~= 4 || size(M, 3) ~= 4
  error('ptt:calibrateChannels:shape', 'M must be [Nt x 4 x 4]');
end
if numel(z) ~= size(M, 1)
  error('ptt:calibrateChannels:depth', 'z has %d rows for %d moment rows', numel(z), size(M, 1));
end

% --- reciprocity: |b/a| and arg a - arg b, summed over the record
mr = z >= recip(1) & z <= recip(2) & isfinite(real(M(:, 3, 3))) & isfinite(real(M(:, 4, 4)));
if nnz(mr) < 10
  error('ptt:calibrateChannels:recip', 'fewer than 10 rows in the reciprocity range [%g %g] m', recip(1), recip(2));
end
p_hv = sum(real(M(mr, 3, 3))); p_vh = sum(real(M(mr, 4, 4)));
x = sum(M(mr, 3, 4));                   % <HV VH*> = a conj(b) |S_x|^2
ba = sqrt(p_vh / max(p_hv, realmin));   % |b| / |a|
dphi = angle(x);                         % arg a - arg b
xcoh = abs(x) / max(sqrt(p_hv * p_vh), realmin);

% --- the firn: |ab| from the co-pol power ratio
mf = z >= firn(1) & z <= firn(2) & isfinite(real(M(:, 1, 1))) & isfinite(real(M(:, 2, 2)));
if nnz(mf) < 10
  error('ptt:calibrateChannels:firn', 'fewer than 10 rows in the firn window [%g %g] m', firn(1), firn(2));
end
firn_ratio = sum(real(M(mf, 2, 2))) / max(sum(real(M(mf, 1, 1))), realmin);
abm = sqrt(firn_ratio);                  % |a b|

% --- the firn: arg(ab) from the co-pol phase at z = 0
% The phase is read per row as the argument of <HH VV*>; the single-axis
% form is fitted through the window on a grid in (mu, rate) with the
% offset phi0 in closed form as the circular mean of the measured-minus-
% model phasors, each row weighted by its coherence. The straight line
% is fitted beside it for the record.
zf = z(mf);
cf = M(mf, 1, 2) ./ max(sqrt(real(M(mf, 1, 1)) .* real(M(mf, 2, 2))), realmin);
w = abs(cf).^2;
[fit_int, fit_dl, fit_mu, rms] = H_firn_phase(cf, zf, w, gpd);
ph = unwrap(angle(cf));
A = [ones(size(zf)), zf];
sw = sqrt(w);                            % weighted LS: residuals weighted by w
coef = (A .* sw) \ (ph .* sw);
line_int = coef(1);
if isempty(phase_ab)
  ab_arg = -fit_int;                     % measured = true - arg(ab), true -> 0 at z = 0
else
  ab_arg = phase_ab;
end
ab_arg = pi - mod(pi - ab_arg, 2*pi);   % onto (-pi, pi], so -pi and pi are one pin

% --- assemble a and b, on the branch in phase with H (see THE SIGN OF a
% AND b): the half-angle split fixes the pair only up to adding pi to both
am = sqrt(abm / ba); bm = sqrt(abm * ba);
arg_a = (ab_arg + dphi) / 2; arg_b = (ab_arg - dphi) / 2;
if cos(arg_a) + cos(arg_b) < 0
  arg_a = arg_a + pi; arg_b = arg_b + pi;
end
arg_a = angle(exp(1i * arg_a)); arg_b = angle(exp(1i * arg_b));
a = am * exp(1i * arg_a); b = bm * exp(1i * arg_b);
g = [1, a * b, a, b];

Mc = M;
for k = 1:4
  for l = 1:4
    Mc(:, k, l) = M(:, k, l) / (g(k) * conj(g(l)));
  end
end

info = struct('a_db', 20*log10(am), 'b_db', 20*log10(bm), ...
  'arg_a_deg', rad2deg(arg_a), 'arg_b_deg', rad2deg(arg_b), ...
  'xpol_phase_deg', rad2deg(dphi), 'xpol_coh', xcoh, ...
  'firn_db', 10*log10(firn_ratio), ...
  'phase_fit', struct('z', zf, 'phase_deg', rad2deg(ph), ...
    'intercept_deg', rad2deg(fit_int), 'dlam', fit_dl, 'mu', fit_mu, ...
    'rms_deg', rad2deg(rms), 'line_intercept_deg', rad2deg(line_int)), ...
  'phase_ab_deg', rad2deg(ab_arg), 'phase_ab_pinned', ~isempty(phase_ab), ...
  'firn_m', firn, 'recip_m', recip, 'g', g);
end

function [phi0, dl, mu_best, rms] = H_firn_phase(cf, zf, w, gpd)
%H_FIRN_PHASE The co-pol phase's value at z = 0 from the single-axis form.
% For each (mu, rate) node the model phase is atan2(mu sin d, A + B cos d)
% with d = rate z; the offset is the weighted circular mean of the
% measured phase less the model, and the node with the smallest weighted
% circular misfit wins. mu = cos 2(psi - theta) spans [-1, 1]; the rate
% spans dlam 0..0.3.
mus = linspace(-1, 1, 41);
dls = linspace(0, 0.3, 61);
best = inf; phi0 = NaN; dl = NaN; mu_best = NaN;
ang = angle(cf);
for mu = mus
  Ak = (1 - mu^2) / 2; Bk = (1 + mu^2) / 2;
  for d = dls
    del = gpd * d * zf;
    pm = atan2(mu * sin(del), Ak + Bk * cos(del));
    r = sum(w .* exp(1i * (ang - pm)));
    p0 = angle(r);
    cost = sum(w .* (1 - cos(ang - pm - p0)));      % circular misfit
    if cost < best, best = cost; phi0 = p0; dl = d; mu_best = mu; end
  end
end
del = gpd * dl * zf;
pm = atan2(mu_best * sin(del), (1 - mu_best^2)/2 + (1 + mu_best^2)/2 * cos(del));
e = angle(exp(1i * (ang - pm - phi0)));
rms = sqrt(sum(w .* e.^2) / max(sum(w), realmin));
end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
