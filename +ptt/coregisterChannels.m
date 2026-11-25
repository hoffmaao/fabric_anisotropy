function [T, info] = coregisterChannels(S, ref, opts)
%COREGISTERCHANNELS Put the polarimetric channels on one range axis: one delay per channel, fitted in the firn.
%
% [T, info] = ptt.coregisterChannels(S, ref, opts)
%
% The channels of a quad-pol record differ by a sub-sample instrument
% delay: uncoregistered, |C_HHVV| is ~0.29 on Ridge A, and an early
% analysis drew a conclusion about antenna geometry from that number. The
% delay is ONE number per channel, so one number per channel is fitted,
% over opts.rows, and the whole channel is shifted by it as an exact
% band-limited shift.
%
% WHY NOT A FIELD OF OFFSETS. The usual coregistration (the OPR toolbox's,
% which built the shipped CSARP_polarimetric products) estimates an offset
% on tiles and lets it vary over the image. In birefringent ice the slow
% and fast echoes arrive tau = delta / (2 pi f) apart (~5 ns at 1500 m on
% Ridge A), so HH and VV, which mix the two with different weights, look
% offset in range by an amount that changes with depth and with the
% axis-to-antenna angle. A tile search measures that as misregistration and
% removes it: on Ridge A the median VV offset it finds is +0.67 samples on
% every N-S frame and -0.24 on the E-W frames, where a true instrument
% delay would be one number (~+0.16). With the band off-centre from the
% baseband reference, each shift also turns the channel phase, and a fabric
% fit reads the result as depth bands in strength (~0.002-0.004) and axis
% (~1.4 deg rms), measured on paired synthetics with a smooth fabric. So
% the delay is fitted where the ice is nearly isotropic, the firn, and
% applied everywhere. (A delay fitted over the WHOLE image lands on neither
% the instrument nor the ice, and was worse on every pair.)
%
% HOW A QUAD-POL SET IS SOLVED. Through its chains, the instrument model
% Rx S Tx: a channel's delay is its receive chain's plus its transmit
% chain's, so relative to hh
%   vv = a + b,  hv = a,  vh = b
% with a the V-receive and b the V-transmit chain's delay. Two pairs pin
% them, each measured where its two channels share their speckle:
%   vv - hh over opts.rows (the firn: isotropic, so hh and vv see the same
%     scatterer, while deeper the birefringence separates them), and
%   vh - hv over opts.rows_x (default every row: reciprocity makes the
%     ice's hv and vh the same at any depth, birefringence or not).
% Any other set is aligned channel by channel against ref over opts.rows.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels (a field named heading is
%         passed through untouched)
%   ref   reference channel name (default 'hh')
%   opts  rows       logical or index rows to estimate from (required: the
%                    firn)
%         rows_x     the rows for the hv-vh pair (default all)
%         search_t   (5) the largest delay searched [samples]
%         oversample (16) the correlation's interpolation factor
%
% Output
%   T     struct of coregistered channels, ref passed through untouched
%   info.row_offset   per-channel delay [samples], positive: the channel
%                     arrives later
%   info.col_offset   0 for every channel (nothing is shifted along track)
%   info.row_med      the same delay, with 0 for ref
%   info.chains       (quad-pol) v_receive and v_transmit, the V chains'
%                     delays against H [samples]

if nargin < 2 || isempty(ref), ref = 'hh'; end
if nargin < 3, opts = struct(); end
mode = H_opt(opts, 'mode', 'constant');
if ~strcmpi(mode, 'constant')
  error('ptt:coregisterChannels:mode', ...
    ['opts.mode is ''%s'': this function fits one delay per channel (''constant'') ' ...
     'and nothing else'], char(mode));
end
rows = H_opt(opts, 'rows', []);
if isempty(rows) || ~any(rows(:))
  error('ptt:coregisterChannels:rows', 'opts.rows is required: the rows to estimate the delay from (the firn)');
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
