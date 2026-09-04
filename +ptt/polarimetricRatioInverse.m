function out = polarimetricRatioInverse(obs, z, opts)
%POLARIMETRICRATIOINVERSE Orientation AND strength from channel power RATIOS.
%
% out = ptt.polarimetricRatioInverse(obs, z, opts)
%
% Inverts the Fujita et al. (2006) model, under a depth-constant fabric
% orientation, for theta0, the eigenvalue difference profile dlam(z) and
% the anisotropic scattering ratio r(z), using the DEPTH-RESOLVED POWER
% RATIOS BETWEEN CHANNELS - principally VV over VH - rather than the
% azimuthal power anomaly of a single channel.
%
% WHY RATIOS, AND WHY THIS PAIR. Two reasons, and the second is the one
% that makes the method usable on an uncalibrated system.
%
% 1. VV/VH CARRIES ORIENTATION AND STRENGTH TOGETHER. With d = theta0 -
%    gamma the angle from the fabric axis to the antenna, r the scattering
%    ratio and delta the accumulated two-way phase difference, the model
%    (see ptt.polarimetricInverse for the reduction) gives
%
%      |s_vv|^2 = r^2 cos^4 d + sin^4 d + 2 r cos^2 d sin^2 d cos(delta)
%      |s_hv|^2 = (1 + r^2 - 2 r cos delta) cos^2 d sin^2 d
%      |s_hh|^2 = cos^4 d + r^2 sin^4 d + 2 r cos^2 d sin^2 d cos(delta)
%
%    so their ratio depends on d through its azimuthal shape and on delta
%    through a depth oscillation, in different functional forms. The two
%    are therefore separable, which is what lets orientation and strength
%    be solved together instead of one being held while the other is
%    fitted. The cross-polarised channel is the right partner because for
%    a single birefringent column its power FACTORISES exactly into an
%    azimuth term sin^2(2d) and a depth term - the co-polarised pair never
%    does (see ptt.quadpolAzimuth).
%
% 2. A CHANNEL GAIN BECOMES AN ADDITIVE CONSTANT. Writing each measured
%    channel as the ice times a transmit and a receive gain, referenced to
%    HH with r_g the V/H receive ratio and t_g the V/H transmit ratio:
%
%      V_hh ~ S_hh          V_hv ~ S_hv t_g
%      V_vh ~ S_vh r_g      V_vv ~ S_vv r_g t_g
%
%    Pairs that SHARE a path cancel one gain outright:
%      VV/VH  share receive V  -> only |t_g|^2 survives
%      HH/HV  share receive H  -> only |t_g|^2 survives
%      HH/VH  share transmit H -> only |r_g|^2 survives
%      HV/VH  share nothing    -> pure instrument, the reciprocity check
%    and in dB a surviving constant gain is a constant OFFSET. So each
%    ratio observable is fitted with its own free additive constant and
%    the calibration drops out of the problem, leaving the SHAPE - in
%    azimuth and in depth - to carry the physics. This is why the method
%    works on a system whose H and V chains differ by 6-29 dB
%    (ptt.equaliseChannels), where every co-polarised-coherence estimator
%    we have is compromised.
%
%    The offsets are constant only where the gains are. They are not
%    across a waveform boundary: a 12.6 dB step was measured at 640-720 m
%    at WAIS Divide, in the transmit chain. Pass opts.offset_edges to give
%    each range block its own offset rather than pretending otherwise.
%
% NO AZIMUTH SYNTHESIS IS USED, AND HEADINGS ARE NOT OPTIONAL. Synthesis
% is the step that would re-mix the channels and destroy the
% offset-cancellation above, so it is avoided: azimuth comes from the
% antenna heading of each column of obs.
%
% A counting argument says one heading suffices - three independent powers
% plus the co-polarised phase against d, r and delta - and that argument
% is misleading. MEASURED in test_polarimetric_ratio: from a single
% heading the strength still comes out (dlam within 0.03) but the
% ORIENTATION lands 18 deg off, because with only one d the fit trades d
% against r(z), which is free at every depth. With four headings the same
% synthetic gives 0.12 deg. So: dlam may be quoted from one heading,
% theta0 may NOT. Supply several headings over the same ice - our grids
% have them - and the axis becomes well determined, because d then varies
% while theta0 does not.
%
% THE 90-DEGREE ALIAS, AND THE ONLY THING THAT BREAKS IT. Power ratios
% alone CANNOT tell (theta0, r) from (theta0 + 90 deg, 1/r). Swapping the
% orientation by a quarter turn exchanges the roles of HH and VV, which
% maps r to 1/r and shifts each ratio by a constant - and a constant is
% exactly what the free gain offsets absorb. Measured on the synthetic in
% test_polarimetric_ratio: the aliased solution fits the four-heading data
% as well as the truth, and adding headings does NOT help, because every
% d = theta0 - gamma shifts together. This is the eigenvector swap Nymand
% warns about, made exact by the offsets.
%
% The CO-POLARISED COHERENCE PHASE breaks it, because under the same swap
% phi changes SIGN, and a sign flip is not something a constant offset can
% absorb. So pass obs.phi wherever the coherence survives and include
% 'phi' in opts.use; its own unknown instrument phase is carried as one
% more free constant. Without it the orientation is determined only to
% within 90 degrees and out.alias_unresolved is set, which callers must
% honour rather than quietly reporting one of the two.
%
% THE CROSS-POLARISED PEDESTAL IS MODELLED, NOT ASSUMED AWAY. This system
% leaks cross-polarised power at about -3.6 dB relative to co-pol, against
% -30 dB for the UWB system this method was developed on. That pedestal
% fills in exactly the nulls the method reads, so it is carried as an
% additive term, ped * (|s_hh|^2 + |s_vv|^2)/2, with ped a free parameter.
% Without it the pedestal is read as a floor on (1 + r^2 - 2r cos delta),
% which biases r high and dlam low.
%
%
% SIGN CONVENTION, PINNED BY TEST. Here d = theta0 - gamma, the SYNTHESIS
% sense, the same one ptt.quadpolAzimuth and ptt.quadpolFabricLS use, in
% which a fabric at azimuth theta0 puts its features at sweep index
% +theta0. ptt.fujitaModel deliberately uses the paper's R S R' sense, in
% which they appear at -theta0. The two agree only with theta negated:
% measured, matching the sign gives correlation 0.967 and median 0.23 dB
% against ptt.fujitaModel, mismatching it gives 0.11 and 1.96 dB. Both are
% asserted in test_polarimetric_inverse, because a silent flip here would
% negate every axis this code reports.
% Inputs
%   obs  .P_hh, .P_vv, .P_hv  [Nz x Ng] measured POWERS (linear, not dB),
%          one column per antenna azimuth. .P_vh optional; if given it is
%          averaged with P_hv AFTER each is offset-fitted, never before.
%        .gamma  [1 x Ng] antenna azimuth of each column (radians,
%          geographic). Scalar allowed for a single heading.
%        .phi    [Nz x Ng] co-pol coherence phase (rad), optional
%   z    [Nz x 1] depth (m)
%   opts .fc (750e6)
%        .n_r (8)      scattering-ratio nodes, piecewise linear in depth
%        .n_d (12)     dlam nodes, piecewise linear in depth
%        .theta0 ([])  initial orientation; [] scans
%        .dlam0 (0.05), .r0 (1), .ped0 (0.3)
%        .eta_r (1e-2), .eta_d (1e-2)  second-difference smoothness
%        .offset_edges ([]) depths splitting the free offsets
%        .use ({'vv_vh','hh_hv'}) any of vv_vh, hh_hv, hh_vh, hh_vv, phi
%        .max_iter (40), .clip_db (-25), .w_phi (2, dB per radian)
%        .n_start (4)  how many of the best grid points to run to
%                 convergence; the loss is multi-modal in BOTH theta0 and
%                 dlam, so one start is not enough
%        .dlam_grid    starting dlam values scanned (default 7 values
%                 spanning 0.01-0.11); pass .dlam0 to pin it
%
% Output
%   .theta0, .dlam_z, .r_z, .ped, .offsets, .loss, .converged, .pred,
%   .rms, .reciprocity_db (measured HV/VH imbalance, a pure instrument
%   diagnostic that should be flat in depth if the model is right)
%   .alias_unresolved  true when 'phi' was not used, meaning theta0 is
%   determined only modulo 90 degrees and r may be 1/r
%
% See also ptt.polarimetricInverse, ptt.quadpolFabricLS, ptt.fujitaModel,
%   ptt.equaliseChannels.

if nargin < 3, opts = struct(); end
z = z(:); Nz = numel(z);
fc = H_opt(opts, 'fc', 750e6);
n_r = H_opt(opts, 'n_r', 8);
n_d = H_opt(opts, 'n_d', 12);
eta_r = H_opt(opts, 'eta_r', 1e-2);
eta_d = H_opt(opts, 'eta_d', 1e-2);
max_iter = H_opt(opts, 'max_iter', 40);
clip_db = H_opt(opts, 'clip_db', -25);
use = H_opt(opts, 'use', {'vv_vh', 'hh_hv'});
if ischar(use), use = {use}; end
% phi is in radians and the ratios in dB, so it needs a scale to sit in
% the same least-squares. MEASURED, not guessed: on the synthetic in
% test_polarimetric_ratio, w_phi = 2 recovers theta0 to 0.12 deg and dlam
% exactly, while 8 and 20 drag the solution to 4.0 and 1.2 deg with the
% pedestal biased low. The phase only has to break the 90-degree alias -
% it is one bit of information - and weighting it to compete with the
% ratios lets its own model error pull the whole fit.
w_phi = H_opt(opts, 'w_phi', 2);

gam = obs.gamma(:).';
Ng = numel(gam);
P = struct('hh', obs.P_hh, 'vv', obs.P_vv, 'hv', obs.P_hv);
if isfield(obs, 'P_vh') && ~isempty(obs.P_vh), P.vh = obs.P_vh; else, P.vh = []; end
for f = {'hh','vv','hv'}
  if size(P.(f{1}), 1) ~= Nz || size(P.(f{1}), 2) ~= Ng
    error('ptt:polarimetricRatioInverse:shape', ...
      'obs.P_%s is %s, expected [%d x %d]', f{1}, mat2str(size(P.(f{1}))), Nz, Ng);
  end
end

% reciprocity diagnostic: HV against VH is pure instrument (r_g/t_g), so
% it should not vary with depth. Reported, never fitted.
recip = [];
if ~isempty(P.vh)
  recip = 10*log10(max(P.vh, realmin) ./ max(P.hv, realmin));
end

% --- offset blocks in depth
edges = H_opt(opts, 'offset_edges', []);
blk = ones(Nz, 1);
if ~isempty(edges)
  for e = edges(:).'
    blk = blk + (z > e);
  end
end
n_blk = max(blk);

C = ptt.constants();
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);
dz = [0; diff(z)];
z_rn = linspace(z(1), z(end), n_r).';
z_dn = linspace(z(1), z(end), n_d).';

% --- pack observed ratios (dB)
phi_obs = H_opt(obs, 'phi', []);
[d_obs, idx, tags, is_ph] = H_pack_obs(P, use, clip_db, phi_obs, w_phi);
if isempty(d_obs)
  error('ptt:polarimetricRatioInverse:empty', 'no finite ratio observations');
end
n_off = numel(use) * n_blk;

% m = [theta0; r_nodes; dlam_nodes; ped; offsets]
th0 = H_opt(opts, 'theta0', []);
m0 = [0; H_opt(opts,'r0',1)*ones(n_r,1); H_opt(opts,'dlam0',0.05)*ones(n_d,1); ...
      H_opt(opts,'ped0',0.3); zeros(n_off,1)];
ip = struct('th', 1, 'r', 1+(1:n_r), 'd', 1+n_r+(1:n_d), ...
  'ped', 2+n_r+n_d, 'off', 2+n_r+n_d+(1:n_off));

fwd = @(mm) H_pack_pred(H_model(mm, gam, z, dz, gpd, z_rn, z_dn, ip, blk, n_blk, use, w_phi), ...
  use, idx, clip_db);

% --- regularisation on the two profiles
Gam = zeros(0, numel(m0));
Gam = [Gam; H_d2(n_r, numel(m0), ip.r) * eta_r];
Gam = [Gam; H_d2(n_d, numel(m0), ip.d) * eta_d];
GtG = Gam.' * Gam;
% RESIDUAL, not a plain difference: the phi entries are ANGLES (scaled into
% the dB least-squares by w_phi), and the instrument phase offset this file
% carries as a free constant puts them near the +-pi branch cut routinely.
% Differencing there scores a true misfit of ~0 as ~2*pi*w_phi and lets
% those entries drag theta0 - the failure ptt.fabricGLS measured (its
% wrapping took chi2/dof from ~2200 to order 1) and fixes the same way.
% out.rms already wrapped; now the objective agrees with the diagnostic.
resid = @(mm) H_resid(d_obs, fwd(mm), is_ph, w_phi);
loss = @(mm) sum(resid(mm).^2) + sum((Gam*mm).^2);

% --- MULTI-START, because this loss is strongly multi-modal and a single
% Gauss-Newton lands wherever it was pointed. Two separate reasons:
%
%   theta0 is periodic, and swapping it by 90 deg exchanges the roles of
%   HH and VV, which the free per-observable offsets can absorb almost
%   exactly - so the eigenvector swap Nymand warns about is not merely a
%   risk here, it is a near-exact alias of the data.
%
%   dlam enters only through cos(delta) with delta ACCUMULATED, so the
%   misfit oscillates in it. Measured on the synthetic in
%   test_polarimetric_ratio: the converged loss varies four-fold across
%   starting dlam between 0.01 and 0.10, and only the start nearest the
%   true shallow value recovers the profile.
%
% So a coarse grid over (theta0, dlam0) is scored, the best n_start
% distinct points are each run to convergence, and the lowest final loss
% wins. Anything less than this reports a local minimum as an answer.
n_start = H_opt(opts, 'n_start', 4);
th_grid = (0:6:174) * pi/180;
if ~isempty(th0), th_grid = th0; end
dl_grid = H_opt(opts, 'dlam_grid', [0.01 0.02 0.03 0.045 0.06 0.08 0.11]);
if isfield(opts, 'dlam0') && ~isempty(opts.dlam0), dl_grid = opts.dlam0; end
cand = []; cl = [];
for th = th_grid
  for d0 = dl_grid
    mt = m0; mt(ip.th) = th; mt(ip.d) = d0;
    cand(end+1, :) = [th d0]; %#ok<AGROW>
    cl(end+1) = loss(mt); %#ok<AGROW>
  end
end
[~, ord] = sort(cl);
n_try = min(n_start, numel(ord));
best_m = []; best_L = inf; best_hist = []; best_conv = false; best_it = 0;
for q = 1:n_try
  m = m0; m(ip.th) = cand(ord(q), 1); m(ip.d) = cand(ord(q), 2);
  [m, L, conv, it] = H_solve(m, resid, GtG, loss, ip, max_iter);
  if L(end) < best_L
    best_L = L(end); best_m = m; best_hist = L; best_conv = conv; best_it = it;
  end
end
m = best_m; L = best_hist; converged = best_conv; it = best_it;
out_starts = n_try;

M = H_model(m, gam, z, dz, gpd, z_rn, z_dn, ip, blk, n_blk, use, w_phi);
out = struct('theta0', mod(m(ip.th), pi), ...
  'r_z', interp1(z_rn, m(ip.r), z, 'linear', 'extrap'), ...
  'dlam_z', max(interp1(z_dn, m(ip.d), z, 'linear', 'extrap'), 0), ...
  'ped', m(ip.ped), 'offsets', reshape(m(ip.off), n_blk, numel(use)), ...
  'loss', L, 'n_iter', it, 'converged', converged, 'pred', M, ...
  'z_rn', z_rn, 'z_dn', z_dn, 'used', {use}, 'gamma', gam, ...
  'reciprocity_db', recip, 'clip_db', clip_db, 'n_starts', out_starts, ...
  'alias_unresolved', ~any(strcmp(use, 'phi')));
out.rms = struct();
for k = 1:numel(use)
  Mk = M.(use{k});
  if ~strcmp(use{k}, 'phi'), Mk = max(Mk, clip_db); end
  R = H_ratio_db(P, use{k}, clip_db, phi_obs, w_phi) - Mk;
  if strcmp(use{k}, 'phi'), R = w_phi * angle(exp(1i*R/w_phi)); end
  out.rms.(matlab.lang.makeValidName(use{k})) = sqrt(mean(R(isfinite(R)).^2));
end
if ~isempty(recip)
  out.reciprocity_spread_db = std(recip(isfinite(recip)));
end
end

% -------------------------------------------------------------------------
function [m, L, converged, it] = H_solve(m, resid, GtG, loss, ip, max_iter)
L = nan(max_iter+1, 1); L(1) = loss(m);
converged = false; it = 0;
for it = 1:max_iter
  r0 = resid(m);
  G = H_jac(resid, m, numel(r0));
  A = G.'*G + GtG;
  b = G.'*r0 - GtG*m;
  dm = A \ b;
  if ~all(isfinite(dm)), break; end
  a = 1; ok = false;
  for h = 1:20
    mt = H_clamp(m + a*dm, ip);
    Lt = loss(mt);
    if Lt < L(it), m = mt; L(it+1) = Lt; ok = true; break; end
    a = a/2;
  end
  if ~ok, converged = true; L(it+1) = L(it); break; end
  if abs(L(it)-L(it+1))/max(L(it+1),realmin) < 1e-5, converged = true; break; end
end
L = L(1:min(it+1, numel(L)));
end

function M = H_model(m, gam, z, dz, gpd, z_rn, z_dn, ip, blk, n_blk, use, w_phi)
th0 = m(ip.th);
r = max(interp1(z_rn, m(ip.r), z, 'linear', 'extrap'), 1e-4);
dl = max(interp1(z_dn, m(ip.d), z, 'linear', 'extrap'), 0);
ped = max(m(ip.ped), 0);
delta = cumsum(dl .* dz) * gpd;                    % [Nz x 1]
d = th0 - gam;                                     % [1 x Ng]
c2 = cos(d).^2; s2 = sin(d).^2;
cd = cos(delta);
% unit-scale channel powers under the Fujita reduction at theta(z)=theta0
Phh = c2.^2 + (r.^2).*s2.^2 + 2*r.*c2.*s2.*cd;
Pvv = (r.^2).*c2.^2 + s2.^2 + 2*r.*c2.*s2.*cd;
Pxx = (1 + r.^2 - 2*r.*cd) .* (c2 .* s2);
% additive antenna-frame cross-pol leakage, as a fraction of co-pol power
Pxx = Pxx + ped * (Phh + Pvv) / 2;
M = struct('Phh', Phh, 'Pvv', Pvv, 'Pxx', Pxx);
off = reshape(m(ip.off), n_blk, numel(use));
for k = 1:numel(use)
  switch use{k}
    case 'vv_vh', R = 10*log10(max(Pvv,realmin)./max(Pxx,realmin));
    case 'hh_hv', R = 10*log10(max(Phh,realmin)./max(Pxx,realmin));
    case 'hh_vh', R = 10*log10(max(Phh,realmin)./max(Pxx,realmin));
    case 'hh_vv', R = 10*log10(max(Phh,realmin)./max(Pvv,realmin));
    case 'phi'
      % Nymand eq 6.10, scaled into the dB least-squares by w_phi. This is
      % the term that changes sign under the 90-deg swap and so is the one
      % that resolves the alias.
      t2 = (s2 ./ max(c2, realmin));            % tan^2 d, [1 x Ng]
      t4 = t2.^2;
      R = w_phi * atan2(r .* sin(delta) .* (1 - t4), ...
                        r .* cos(delta) .* (1 + t4) + t2 .* (1 + r.^2));
    otherwise, error('unknown observable %s', use{k});
  end
  M.(use{k}) = R + off(blk, k);        % free additive gain offset per block
end
end

function [d, idx, tags, is_ph] = H_pack_obs(P, use, clip_db, phi, w_phi)
d = []; idx = struct(); tags = {}; is_ph = [];
for k = 1:numel(use)
  R = H_ratio_db(P, use{k}, clip_db, phi, w_phi);
  g = isfinite(R);
  idx.(matlab.lang.makeValidName(use{k})) = g;
  d = [d; R(g)]; %#ok<AGROW>
  is_ph = [is_ph; repmat(strcmp(use{k}, 'phi'), nnz(g), 1)]; %#ok<AGROW>
  tags{end+1} = use{k}; %#ok<AGROW>
end
is_ph = logical(is_ph);
end

function R = H_ratio_db(P, name, clip_db, phi, w_phi)
if strcmp(name, 'phi')
  if nargin < 4 || isempty(phi)
    error('ptt:polarimetricRatioInverse:phi', ...
      "opts.use lists 'phi' but obs.phi was not supplied");
  end
  R = w_phi * phi;
  return;
end
switch name
  case 'vv_vh'
    den = P.vh; if isempty(den), den = P.hv; end
    R = 10*log10(max(P.vv, realmin) ./ max(den, realmin));
  case 'hh_hv', R = 10*log10(max(P.hh, realmin) ./ max(P.hv, realmin));
  case 'hh_vh'
    den = P.vh; if isempty(den), den = P.hv; end
    R = 10*log10(max(P.hh, realmin) ./ max(den, realmin));
  case 'hh_vv', R = 10*log10(max(P.hh, realmin) ./ max(P.vv, realmin));
  otherwise, error('unknown observable %s', name);
end
R(~isfinite(R)) = NaN;
R = max(R, clip_db);
end

function v = H_pack_pred(M, use, idx, clip_db)
% The floor is applied to the MODEL as it is to the data (H_ratio_db).
% Flooring only the data leaves the model free to run to 10log10(realmin)
% at its own power nulls, generating hundreds of dB of residual at exactly
% the cells the clip exists to neutralise, and the fit then trades
% theta0/dlam to move a null whose depth it cannot know. 'phi' is a phase,
% not a ratio, and is never floored.
v = [];
for k = 1:numel(use)
  R = M.(use{k});
  if ~strcmp(use{k}, 'phi'), R = max(R, clip_db); end
  g = idx.(matlab.lang.makeValidName(use{k}));
  v = [v; R(g)]; %#ok<AGROW>
end
end

function D = H_d2(n, ncol, cols)
D = zeros(max(n-2,0), ncol);
for k = 1:n-2, D(k, cols(k:k+2)) = [1 -2 1]; end
end

function m = H_clamp(m, ip)
m(ip.r) = max(m(ip.r), 1e-3);
m(ip.d) = max(m(ip.d), 0);
m(ip.ped) = min(max(m(ip.ped), 0), 5);
end

function G = H_jac(resid, m, nd)
% Jacobian of the FORWARD, differenced through the residual so the phase
% wrap is inside it: d(resid)/dm = -d(fwd)/dm away from the branch cut, and
% at the cut the wrapped difference is the one the objective actually uses.
np = numel(m); G = zeros(nd, np); r0 = resid(m);
for j = 1:np
  h = max(1e-7, 1e-4*abs(m(j)));
  mp = m; mp(j) = mp(j) + h;
  G(:, j) = -(resid(mp) - r0) / h;
end
G(~isfinite(G)) = 0;
end

function r = H_resid(d, g, is_ph, w_phi)
r = d - g;
r(is_ph) = w_phi * angle(exp(1i * r(is_ph) / w_phi));
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
