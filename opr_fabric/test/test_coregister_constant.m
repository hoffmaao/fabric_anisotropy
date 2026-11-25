%TEST_COREGISTER_CONSTANT One delay per channel, from the firn, applied everywhere.
%
% ptt.coregisterChannels(S, ref, struct('mode', 'constant', 'rows', R))
% fits ONE range delay per channel over the rows R and shifts the whole
% channel by it. The tile mode follows any delay that changes with depth -
% in birefringent ice that includes the HH-VV lag of the slow and fast
% echoes, which is signal - so the constant mode is estimated where the ice
% is nearly isotropic (the firn) and leaves the deeper lag alone.
%
% Synthetic: band-limited complex speckle (Hann spectrum, 160 MHz at
% 300 MHz sampling, as the accumulation radar), 4000 rows by 300 traces,
% under one layered power envelope that every channel shares (the
% stratigraphy). The instrument is a V-receive chain 0.12 samples and a
% V-transmit chain 0.25 samples later than H, so VV is HH delayed by 0.37,
% HV by 0.12 and VH by 0.25 (HV and VH the same cross-pol speckle, as
% reciprocity makes them); below row 800 VV is delayed by a further 0 to
% 1.2 samples, growing with depth (the ice's birefringent lag).
%
%   1. ESTIMATE. With rows 1-600 as the firn, the three delays come back
%      within 0.02 samples, and the deeper lag does not pull them.
%   2. APPLY. After the shift VV matches HH in rows 1-600 (coherence
%      > 0.99) while the deep lag is still there (coherence well below).
%   3. BLANKED SAMPLES stay blanked where they land and do not spread.
%   4. ROWS ARE REQUIRED: without opts.rows the call errors.
%
% Fast tier (seconds). Run: matlab -batch "run('opr_fabric/test/test_coregister_constant.m')"
clear;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));
fails = 0;
rs = RandStream('mt19937ar', 'Seed', 7);
Nt = 4000; Nx = 300; N = 8192;
f = [0:N / 2 - 1, -N / 2:-1].' / N;                       % cycles per sample
W = zeros(N, 1); band = abs(f) <= 0.5 * 160 / 300;       % 160 MHz of 300
W(band) = hann(nnz(band));
spk = @() ifft(W .* (randn(rs, N, Nx) + 1i * randn(rs, N, Nx)));
delay = @(X, d) ifft(fft(X) .* exp(-1i * 2 * pi * f .* d));  % later by d samples (d may vary by row: done per segment below)
% layers: a log-normal reflectivity varying over ~10 samples, common to all
% channels and traces
env = exp(0.8 * conv(randn(rs, N, 1), hann(21) / sum(hann(21)) * sqrt(21), 'same'));
A = spk() .* env; B = 0.5 * spk() .* env;
hh = A(1:Nt, :);
% VV: the constant 0.37 everywhere, plus a lag growing from 0 at row 800 to
% 1.2 samples at the bottom, built segment by segment (each 200 rows held)
vv = zeros(Nt, Nx);
extra = max(0, ((1:Nt).' - 800) / (Nt - 800)) * 1.2;
for r0 = 1:200:Nt
  r = r0:min(Nt, r0 + 199);
  Y = delay(A, 0.37 + mean(extra(r)));
  vv(r, :) = Y(r, :);
end
Y = delay(B, 0.12); hv = Y(1:Nt, :);
Y = delay(B, 0.25); vh = Y(1:Nt, :);
S = struct('hh', hh, 'vv', vv, 'hv', hv, 'vh', vh);
rows = (1:Nt).' <= 600;

%% 1. estimate
[T, info] = ptt.coregisterChannels(S, 'hh', struct('mode', 'constant', 'rows', rows));
e = [info.row_med.vv - 0.37, info.row_med.hv - 0.12, info.row_med.vh - 0.25];
ok = all(abs(e) < 0.02) && info.row_med.hh == 0;
fprintf('1. delays: VV %+.3f (true +0.370), HV %+.3f (+0.120), VH %+.3f (+0.250); chains V-rx %+.3f, V-tx %+.3f: %s\n', ...
  info.row_med.vv, info.row_med.hv, info.row_med.vh, info.chains.v_receive, info.chains.v_transmit, H_tick(ok));
fails = fails + ~ok;

%% 2. apply
coh = @(a, b, r) abs(sum(sum(a(r, :) .* conj(b(r, :))))) / sqrt(sum(sum(abs(a(r, :)).^2)) * sum(sum(abs(b(r, :)).^2)));
c_top = coh(T.hh, T.vv, 50:550); c_deep = coh(T.hh, T.vv, 3500:3950); c_x = coh(T.hv, T.vh, 50:3950);
ok = c_top > 0.99 && c_deep < 0.95 && c_x > 0.99;
fprintf('2. after the shift: HH-VV coherence %.4f in the firn rows, %.3f deep (the lag is left); HV-VH %.4f: %s\n', ...
  c_top, c_deep, c_x, H_tick(ok));
fails = fails + ~ok;

%% 3. blanked samples
Sb = S; Sb.vv(1000:1010, 5) = NaN;
Tb = ptt.coregisterChannels(Sb, 'hh', struct('mode', 'constant', 'rows', rows));
nb = nnz(isnan(Tb.vv)); hit = find(isnan(Tb.vv(:, 5)));
ok = nb == 11 && all(abs(hit - (1000:1010).') <= 1) && all(isfinite(Tb.vv(:, [1:4, 6:end])), 'all');
fprintf('3. 11 blanked samples come back as %d, at rows %d-%d of trace 5 only: %s\n', nb, min(hit), max(hit), H_tick(ok));
fails = fails + ~ok;

%% 4. rows are required
try
  ptt.coregisterChannels(S, 'hh', struct('mode', 'constant'));
  ok = false;
catch err
  ok = strcmp(err.identifier, 'ptt:coregisterChannels:rows');
end
fprintf('4. no rows -> error ptt:coregisterChannels:rows: %s\n', H_tick(ok));
fails = fails + ~ok;

if fails > 0
  error('test_coregister_constant: %d check(s) failed', fails);
end
fprintf('test_coregister_constant: all checks passed\n');

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
