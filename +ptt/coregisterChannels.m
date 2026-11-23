function [T, info] = coregisterChannels(S, ref, opts)
%COREGISTERCHANNELS Align polarimetric channels using the OPR coregistration.
%
% [T, info] = ptt.coregisterChannels(S, ref, opts)
%
% Thin wrapper around the CReSIS/OPR toolbox function `coregistration`
% (opr/matlab/processing/coregistration.m), applied to every channel of S
% against `ref` (default 'hh').
%
% CALLS THE TOOLBOX RATHER THAN REIMPLEMENTING IT. The production
% polarimetric products are built by polarimetric_task.m calling exactly
% this function, so using it here means the quad-pol path is aligned by the
% same code, with the same defaults, as the CSARP_polarimetric products
% every other result in this repo is derived from. A hand-rolled
% cross-spectrum fit was tried first and was worse on both counts: it
% reached |C_HHVV| 0.371 against the production 0.495, and it DEGRADED the
% already-aligned cross-polarized pair from 0.985 to 0.882, because a
% single global delay is the wrong model for an offset that varies over the
% image. (Much of that variation is the ice's birefringent lag, see below,
% so a delay fitted over the WHOLE image lands on neither the instrument nor
% the ice; mode 'constant' fits it in the firn and through the receive and
% transmit chains instead.) The toolbox version estimates offsets on overlapping tiles in both
% range and azimuth, oversamples the correlation, median-filters the offset
% field and resamples by sinc interpolation.
%
% WHY IT MATTERS HERE. Reading the standardphase channels raw and skipping
% this step gives |C_HHVV| ~ 0.29 on Ridge A, against 0.495 for the
% coregistered product - and an earlier analysis in this repo drew a
% conclusion about antenna geometry from that uncorrected number. The
% channels are offset by a sub-sample instrument delay, not by looking at
% different ice.
%
% BUT THE TILES ALSO FOLLOW THE BIREFRINGENCE (opts.mode, 7 Nov 2026). In
% birefringent ice the slow and fast echoes arrive tau = delta / (2 pi f)
% apart (~5 ns at 1500 m on Ridge A), so HH and VV, which mix the two with
% different weights, look offset in range by an amount that changes with
% depth and with the axis-to-antenna angle. The tile search measures that as
% misregistration and removes it: on Ridge A the median VV offset is +0.67
% samples on every N-S frame and -0.24 on the E-W frames, where a true
% instrument delay would be one number (~+0.16). With the band off-centre
% from the baseband reference, each shift also turns the channel phase, and
% the fabric fits read the result as depth bands in strength (~0.002-0.004)
% and axis (~1.4 deg rms), measured on paired synthetics with a smooth
% fabric. mode 'constant' fits ONE delay per channel over opts.rows - the
% firn, where the ice is nearly isotropic - and applies it everywhere.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels
%   ref   reference channel name (default 'hh')
%   opts  mode ('tiles' default: the OPR tile field below; 'constant': one
%           delay per channel, from rows opts.rows, see above)
%         for 'tiles', passed to `coregistration`; defaults below match
%           polarimetric.m: Tt 51, Tx 101, overlap_t 25, overlap_x 25,
%           search_t 5, search_x 5, one_dim_search_en true
%         for 'constant': rows (logical or index rows to estimate from,
%           required: the firn), rows_x (the rows for the hv-vh pair,
%           default all), search_t (5: the largest delay searched,
%           samples), oversample (16: the correlation's interpolation
%           factor); see H_constant for how a quad-pol set is solved
%         impl ('toolbox' default: the OPR function; 'ptt': the package
%           copy ptt.coregistration, identical unless its own options are
%           set) and, for 'ptt' only, peak_clamp_en (false) and
%           peak_contrast_min (0) - see ptt.coregistration for what they
%           fix and for the margin measurement that motivated the copy
%
% Output
%   T     struct of coregistered channels, ref passed through untouched
%   info.row_offset / .col_offset  per-channel offset fields ('constant':
%                                  the one delay, and 0)
%   info.row_med                   median row offset per channel [samples]
%                                  (positive: the channel arrives later)
%   info.chains ('constant', quad-pol) v_receive and v_transmit, the V
%                                  chains' delays against H [samples]

if nargin < 2 || isempty(ref), ref = 'hh'; end
if nargin < 3, opts = struct(); end
if strcmpi(H_opt(opts, 'mode', 'tiles'), 'constant')
  [T, info] = H_constant(S, ref, opts);
  return
end
impl = H_opt(opts, 'impl', 'toolbox');
if strcmpi(impl, 'ptt')
  coreg_fn = @ptt.coregistration;
  extra = {'peak_clamp_en', logical(H_opt(opts, 'peak_clamp_en', false)), ...
    'peak_contrast_min', H_opt(opts, 'peak_contrast_min', 0)};
elseif strcmpi(impl, 'toolbox')
  coreg_fn = @coregistration;
  extra = {};
  if exist('coregistration', 'file') ~= 2
    error('ptt:coregisterChannels:toolbox', ...
      ['the OPR toolbox function `coregistration` is not on the path; add ' ...
       'opr/matlab/processing (the startup script does this on the ' ...
       'server) - this wrapper deliberately does not carry its own ' ...
       'implementation']);
  end
else
  error('ptt:coregisterChannels:impl', ...
    ['opts.impl is ''%s''; it must be ''toolbox'' (the OPR function) or ' ...
     '''ptt'' (the package copy ptt.coregistration). Misspelling it must ' ...
     'not silently select the other one, since the two differ in ' ...
     'exactly the peak clamp this option exists to reach.'], impl);
end
if ~isfield(S, ref)
  error('ptt:coregisterChannels:ref', 'reference channel %s not in S', ref);
end

% polarimetric.m's defaults, so this matches how the products were built
Tt = H_opt(opts, 'Tt', 51);
Tx = H_opt(opts, 'Tx', 101);
ov_t = H_opt(opts, 'overlap_t', 25);
ov_x = H_opt(opts, 'overlap_x', 25);
se_t = H_opt(opts, 'search_t', 5);
se_x = H_opt(opts, 'search_x', 5);
one_d = H_opt(opts, 'one_dim_search_en', true);

fn = fieldnames(S);
T = struct();
info = struct('row_offset', struct(), 'col_offset', struct(), ...
  'row_med', struct());
for k = 1:numel(fn)
  ch = fn{k};
  if strcmp(ch, ref)
    T.(ch) = S.(ch);
    info.row_med.(ch) = 0;
    continue;
  end
  [reg, row_off, col_off] = coreg_fn(S.(ref), S.(ch), ...
    'Tt', Tt, 'Tx', Tx, 'overlap_t', ov_t, 'overlap_x', ov_x, ...
    'search_t', se_t, 'search_x', se_x, 'one_dim_search_en', one_d, ...
    'debug_print_en', false, extra{:});
  T.(ch) = reg;
  info.row_offset.(ch) = row_off;
  info.col_offset.(ch) = col_off;
  info.row_med.(ch) = median(row_off(isfinite(row_off)));
end

end

function [T, info] = H_constant(S, ref, opts)
% One delay per channel, applied to the whole record as an exact
% band-limited shift.
%
% A quad-pol set against hh is solved through its chains, the instrument
% model the fabric fits use (Rx S Tx): a channel's delay is its receive
% chain's plus its transmit chain's, so relative to hh
%   vv = a + b,  hv = a,  vh = b
% with a the V-receive and b the V-transmit chain's delay. Two pairs pin
% them, each measured where its two channels share their speckle:
%   vv - hh over opts.rows (the firn: isotropic, so hh and vv see the same
%     scatterer, while deeper the birefringence separates them), and
%   vh - hv over opts.rows_x (default every row: reciprocity makes the
%     ice's hv and vh the same at any depth, birefringence or not).
% Any other set is aligned channel by channel against ref over opts.rows.
rows = H_opt(opts, 'rows', []);
if isempty(rows) || ~any(rows(:))
  error('ptt:coregisterChannels:rows', 'mode ''constant'' needs opts.rows, the rows to estimate the delay from');
end
maxlag = H_opt(opts, 'search_t', 5);
ovs = H_opt(opts, 'oversample', 16);
fn = setdiff(fieldnames(S), {'heading'}, 'stable');
T = S;
info = struct('row_offset', struct(), 'col_offset', struct(), 'row_med', struct());
d = struct();
if strcmp(ref, 'hh') && all(ismember({'hh', 'vv', 'hv', 'vh'}, fn))
  rows_x = H_opt(opts, 'rows_x', true(size(S.hh, 1), 1));
  d_co = H_lag(S.hh(rows, :), S.vv(rows, :), maxlag, ovs);       % a + b
  d_x = H_lag(S.hv(rows_x, :), S.vh(rows_x, :), maxlag, ovs);    % b - a
  d.vv = d_co; d.hv = (d_co - d_x) / 2; d.vh = (d_co + d_x) / 2;
  info.chains = struct('v_receive', d.hv, 'v_transmit', d.vh);
else
  for k = 1:numel(fn)
    if ~strcmp(fn{k}, ref)
      d.(fn{k}) = H_lag(S.(ref)(rows, :), S.(fn{k})(rows, :), maxlag, ovs);
    end
  end
end
for k = 1:numel(fn)
  ch = fn{k};
  if strcmp(ch, ref)
    info.row_med.(ch) = 0;
    continue
  end
  T.(ch) = H_shift(S.(ch), d.(ch));
  info.row_offset.(ch) = d.(ch);
  info.col_offset.(ch) = 0;
  info.row_med.(ch) = d.(ch);
end
end

function lag = H_lag(a, b, maxlag, ovs)
% Delay of b against a in samples (positive: b later): the peak of the
% magnitude of the complex cross-correlation along range, summed over
% traces. The complex correlation is band-limited like the signal, so it
% interpolates exactly (oversampled by ovs, then a parabola) - unlike a
% correlation of |x|^2, whose doubled band aliases at this radar's sampling.
a(~isfinite(a)) = 0; b(~isfinite(b)) = 0;
N = 2 ^ nextpow2(2 * size(a, 1));
c = ifft(sum(conj(fft(a, N)) .* fft(b, N), 2));             % c(k) = sum_n conj(a(n)) b(n + k)
c = fftshift(c);                                             % lag 0 at N/2 + 1
co = abs(interpft(c, ovs * N));                              % periodic, so exact interpolation
lags = ((0:ovs * N - 1).' - ovs * N / 2) / ovs;
in = abs(lags) <= maxlag;
co(~in) = -Inf;
[~, i] = max(co);
if i > 1 && i < numel(co) && all(isfinite(co(i - 1:i + 1)))
  y = co(i - 1:i + 1);
  dd = 0.5 * (y(1) - y(3)) / (y(1) - 2 * y(2) + y(3));
else
  dd = 0;
end
lag = lags(i) + dd / ovs;
end

function y = H_shift(x, lag)
% x advanced by lag samples, y(n) = x(n + lag), by a phase ramp on the
% zero-padded spectrum along range. Blanked (NaN) samples are filled with
% zero for the transform and blanked again where they land.
bad = ~isfinite(x);
x(bad) = 0;
Nt = size(x, 1);
N = 2 ^ nextpow2(Nt + 2 * ceil(abs(lag)) + 64);
kf = [0:N / 2 - 1, -N / 2:-1].';
y = ifft(fft(x, N) .* exp(1i * 2 * pi * kf * lag / N));
y = y(1:Nt, :);
if any(bad(:))
  sh = round(lag);
  badm = false(size(bad));
  src = (1:Nt) + sh;
  okr = src >= 1 & src <= Nt;
  badm(okr, :) = bad(src(okr), :);
  y(badm) = NaN;
end
end

function v = H_opt(o, fld, d)
if isstruct(o) && isfield(o, fld) && ~isempty(o.(fld)), v = o.(fld);
else, v = d; end
end
