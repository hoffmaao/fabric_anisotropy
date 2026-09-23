function [Me, G] = equaliseChannels(M, z, opts)
%EQUALISECHANNELS Range-dependent per-channel gain equalisation of moments.
%
% [Me, G] = ptt.equaliseChannels(M, z, opts)
%
% Rescales the ANTENNA-frame moment matrix M [Nt x 4 x 4] (hh, vv, hv, vh)
% so that the four channels sit on one gain scale at every range, from
% two physical constraints and nothing else:
%
%   RECIPROCITY. The cross-polarised pair must satisfy |HV|^2 = |VH|^2 at
%   every range for a reciprocal medium and antennas. The ratio measured
%   on the data is therefore the gain ratio of the two chains, range bin
%   by range bin, and HV is scaled onto VH (VH and HH are kept as the
%   reference, the H-receive chain).
%
%   THE FIRN AND THE RECEIVER'S RANGE DEPENDENCE. The co-polarised pair
%   has no reciprocity, so VV's offset is anchored in the firn (opts.firn_m,
%   default 160-500 m), where the fabric is weak and VV/HH should be 0 dB
%   up to reflection anisotropy; its RANGE dependence is taken from the
%   cross-pol ratio's, because VV and HV share the V-receive chain and
%   step together at the waveform/gain seam (measured on WAIS Divide:
%   both drop 12.6-12.8 dB relative to HH/VH between 640 and 720 m).
%
% MEASURED (2 Sep 2026, ptt block moments): Ridge A (2024_Antarctica_Ground2)
% is within 3 dB on both ratios and needs nothing; WAIS Divide
% (2023_Antarctica_Ground) has VH/HV = +17 dB above the seam and +29 dB
% below it, VV/HH +22 dB at the surface, -3 dB in the firn, -17 dB deep;
% EastGRIP (2024_Greenland_Ground2) VH/HV +16 -> +12 dB smoothly and
% VV/HH swinging +12 .. -9 dB. Unequalised, WAIS's synthesized co-pol
% power pattern peaks at the H antenna and nulls at the V antenna at
% every depth (modulation 0.98) and every azimuth-synthesising estimator
% reads that as a fabric axis on the antenna; equalised, the modulation
% falls to 0.3-0.6, the extremes move off the antennas and wander with
% depth as fabric nodes do, and the coherence LS held axis lands at
% 113-116 deg with 9 deg spread and dlam 0.11.
%
% This is a moment-level stand-in for the planned raw-channel chan_equal:
% M_kl scales by g_k g_l, exactly what scaling the channels before the
% moments would do, so the same table applies at either level.
%
% Inputs
%   M     [Nt x 4 x 4] antenna-frame moments (ptt.quadpolMoments)
%   z     [Nt x 1] depth (m)
%   opts  firn_m ([160 500]), smooth_m (50, median filter length on the
%         ratio), copol_mode ('firn' default | 'none' leaves VV alone |
%         'reciprocity' applies the full cross-pol gain to VV too, which
%         over-corrects by the co-pol/cross-pol offset - kept for tests)
%
% Outputs
%   Me    equalised moments, same shape
%   G     struct: g [Nt x 4] amplitude gains applied per channel,
%         vhhv_db [Nt x 1] the measured VH/HV before, vvhh_db [Nt x 1]
%         VV/HH before, firn_db the VV anchor applied, seam_db the
%         range-dependence span of the V-receive gain over the record
%
% See also ptt.quadpolMoments, ptt.quadpolFabricPower.

if nargin < 3, opts = struct(); end
firn = H_opt(opts, 'firn_m', [160 500]);
smooth_m = H_opt(opts, 'smooth_m', 50);
mode = H_opt(opts, 'copol_mode', 'firn');
z = z(:);
Nt = numel(z);
dz = median(abs(diff(z)));
nr = max(3, 2*floor(round(smooth_m / max(dz, eps)) / 2) + 1);
kr = ones(nr, 1) / nr;
p = zeros(Nt, 4);
for k = 1:4
  p(:, k) = conv(max(real(M(:, k, k)), 0), kr, 'same');
end
vhhv = p(:, 4) ./ max(p(:, 3), realmin);
vvhh = p(:, 2) ./ max(p(:, 1), realmin);
gHV = medfilt1(sqrt(vhhv), nr);
gHV(~isfinite(gHV) | gHV <= 0) = 1;
inf_ = z >= firn(1) & z <= firn(2);
if nnz(inf_) < 10
  error('ptt:equaliseChannels:firn', 'fewer than 10 rows in the firn window [%g %g] m', firn);
end
switch lower(mode)
  case 'none'
    gVV = ones(Nt, 1);
    g0 = 1;
  case 'reciprocity'
    gVV = gHV;
    g0 = NaN;
  otherwise
    g0 = sqrt(median(p(inf_, 1) ./ max(p(inf_, 2), realmin)));   % VV onto HH in the firn
    step = gHV / median(gHV(inf_));                              % the V receiver's range dependence, unity in the firn
    gVV = g0 * step;
end
g = [ones(Nt, 1), gVV, gHV, ones(Nt, 1)];
Me = M;
for k = 1:4
  for l = 1:4
    Me(:, k, l) = M(:, k, l) .* g(:, k) .* g(:, l);
  end
end
G = struct('g', g, 'vhhv_db', 10*log10(vhhv), 'vvhh_db', 10*log10(vvhh), ...
  'firn_db', 20*log10(g0), 'seam_db', 20*log10(max(gHV) / max(min(gHV), realmin)), ...
  'firn_m', firn, 'mode', mode);
end

function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
