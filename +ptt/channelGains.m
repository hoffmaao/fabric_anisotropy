function [g, info] = channelGains(S, z, opts)
%CHANNELGAINS Complex V-chain gains from reciprocity and the firn, to put four channels on one scale.
%
% [g, info] = ptt.channelGains(S, z, opts)
%
% With H as the reference chain, the measured channels of a reciprocal
% medium are hh = s_hh, vv = gr gt s_vv, hv = gr s_x, vh = gt s_x, where
% gr and gt are the complex gains of the V receive and V transmit chains
% relative to H. Two quantities on the data fix them:
%
%   RECIPROCITY, at every depth: <hv conj(vh)> = gr conj(gt) <|s_x|^2>
%   and |hv|^2 / |vh|^2 = |gr / gt|^2, so gr / gt follows wherever there
%   is cross-polarised power (on this radar at every depth, coherence
%   above 0.98).
%   THE FIRN, nearly isotropic: s_vv = s_hh there, so <vv conj(hh)> /
%   <|hh|^2> = gr gt, phase included. The firn window is kept shallow
%   (opts.firn_m, 30-90 m) because the window's own birefringent phase
%   goes into this estimate; whatever is left there becomes the
%   travel-time reference, which is what the shallow estimators reference
%   to anyway.
%
% gr = sqrt(gr gt * gr / gt) up to a sign, and gt = gr gt / gr. The sign
% flips hv and vh together, which is the mirror of every synthesised
% azimuth about the antenna frame; it is taken with Re gr > 0, the branch
% the Ridge A heading test selected for this instrument.
%
% Why this matters: uncalibrated, the channels synthesised at any
% azimuth are impure - a fraction (1 - gr)(1 - gt) cos^2 psi sin^2 psi
% of the other wave leaks into each - and the eigen-pair phase then
% advances unevenly through each fringe, a strength bias that alternates
% in sign with the node ladder (measured 20-40 % on a synthetic with the
% gains this system shows and an axis 35 deg from the antennas).
%
% Inputs
%   S     struct hh, vv, hv, vh [Nt x Nx] complex, antenna frame
%   z     [Nt x 1] depth [m]
%   opts  firn_m ([30 90]); x_rows (logical [Nt x 1], rows for the
%         reciprocity estimate; default those with cross-pol coherence
%         above coh_min); coh_min (0.8)
%
% Outputs
%   g     struct gr, gt (complex), G = gr gt, R = gr / gt
%   info  coh_firn (co-pol coherence in the firn window), coh_x (cross-pol
%         coherence over the rows used), n_x (rows used), copol_phase_deg,
%         copol_db, xpol_phase_deg (HV minus VH), xpol_db
if nargin < 3, opts = struct(); end
firn = H_opt(opts, 'firn_m', [30 90]);
coh_min = H_opt(opts, 'coh_min', 0.8);
z = double(z(:));
hh = double(S.hh); vv = double(S.vv); hv = double(S.hv); vh = double(S.vh);
% a sample any channel lacks leaves every sum below; zero it in all four, so
% each estimate sums the same samples of both its channels
bad = ~(isfinite(hh) & isfinite(vv) & isfinite(hv) & isfinite(vh));
if any(bad, 'all')
  hh(bad) = 0; vv(bad) = 0; hv(bad) = 0; vh(bad) = 0;
end
rf = z >= firn(1) & z <= firn(2);
if nnz(rf) < 5
  error('ptt:channelGains:firn', 'fewer than 5 rows in the firn window [%g %g] m', firn);
end
c_vh = sum(vv(rf, :) .* conj(hh(rf, :)), 'all');
p_h = sum(abs(hh(rf, :)) .^ 2, 'all'); p_v = sum(abs(vv(rf, :)) .^ 2, 'all');
G = sqrt(p_v / p_h) * exp(1i * angle(c_vh));
coh_firn = abs(c_vh) / sqrt(p_h * p_v);
% reciprocity: rows where the cross-pol pair is coherent
if isfield(opts, 'x_rows') && ~isempty(opts.x_rows)
  rx = logical(opts.x_rows(:));
else
  cx = abs(sum(hv .* conj(vh), 2)) ./ sqrt(sum(abs(hv) .^ 2, 2) .* sum(abs(vh) .^ 2, 2));
  rx = cx >= coh_min;
  if nnz(rx) < 5, rx = isfinite(cx); end
end
c_x = sum(hv(rx, :) .* conj(vh(rx, :)), 'all');
p_hv = sum(abs(hv(rx, :)) .^ 2, 'all'); p_vh = sum(abs(vh(rx, :)) .^ 2, 'all');
R = sqrt(p_hv / p_vh) * exp(1i * angle(c_x));
coh_x = abs(c_x) / sqrt(p_hv * p_vh);
gr = sqrt(G * R);
if real(gr) < 0, gr = -gr; end
gt = G / gr;
g = struct('gr', gr, 'gt', gt, 'G', G, 'R', R);
info = struct('coh_firn', coh_firn, 'coh_x', coh_x, 'n_x', nnz(rx), ...
  'copol_phase_deg', rad2deg(angle(G)), 'copol_db', 20 * log10(abs(G)), ...
  'xpol_phase_deg', rad2deg(angle(R)), 'xpol_db', 20 * log10(abs(R)));
end

function v = H_opt(o, name, default)
if isfield(o, name) && ~isempty(o.(name)), v = o.(name); else, v = default; end
end
