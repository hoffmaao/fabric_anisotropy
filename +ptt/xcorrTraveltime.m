function out = xcorrTraveltime(a, b, dt, opts)
%XCORRTRAVELTIME Travel-time difference between two polarisations by windowed cross-correlation.
%
% out = ptt.xcorrTraveltime(a, b, dt, opts)
%
% Step 2's measurement in the Nymand two-step scheme (the poster's Figure
% 2, after Zeising et al. 2023): in successive fast-time windows the two
% complex signals are cross-correlated over a fine grid of lags, and the
% lag of best correlation is the two-way travel-time difference
%
%   dtau = t_a - t_b   [ns]   (positive: a arrives LATER, a is the slow axis)
%
% between them at that depth. Fed to traveltimeFabricML, with a and b the
% eigen-polarisations of the column, dtau(z) gives dlam(z) by linear
% maximum likelihood.
%
% HOW THE LAG IS FOUND. The signals are complex baseband, so a delay tau
% shows up twice: the ENVELOPE of the cross-correlation peaks at tau with
% a width of about 1/B (6 ns on the accumulation radar), and the carrier
% PHASE at the peak is 2*pi*fc*tau modulo one fringe (1.3 ns). The
% envelope alone is too blunt and the phase alone is ambiguous, so the
% two are combined: the phase gives the fine delay, the fringe integer is
% carried down the column by CONTINUITY from the surface, where the delay
% is zero, and the envelope is the independent check - a window whose
% envelope peak sits more than half a fringe from the tracked delay is
% flagged (out.env_disagree) rather than silently accepted. The tracked
% delay changes by ~0.1 ns per 25 m window on Ridge A against a 1.3 ns
% fringe, so continuity is well determined; a genuine slip would show as
% a run of env_disagree flags below one depth.
%
% WHAT THE SIGNALS MUST BE. a and b must be the SAME scattering seen
% through the two eigen-polarisations - synthesise them with
% ptt.rotatePolarization at the orientation step 1 returns, and give it
% channels whose coregistration did not already remove the birefringent
% delay (one constant delay per channel from the firn, not a per-tile
% offset field that tracks the lag). A pair rotated by dth from the
% eigenframe does NOT read cos(2*dth) of the delay (measured,
% test_xcorr_traveltime verdict 5): it stays near the full delay out to
% ~30 deg and drops to zero at 45 deg, with no loss of coherence.
%
% Inputs
%   a, b   complex [Nt x Nx]: fast time down the columns, traces across.
%          All traces in the array are pooled into one estimate per
%          window - pass one along-track block at a time.
%   dt     fast-time sample interval [ns]
%   opts   (all optional)
%     .fc        (750e6)  the frequency the stored baseband was mixed down
%                         by [Hz] - a delay tau leaves exp(-+2 pi i fc tau)
%                         on a baseband echo, and that phase is the fine
%                         delay. This is the LO, not the band's centroid:
%                         the band sitting off-centre (about -13 MHz on
%                         the accumulation radar) is part of the baseband
%                         spectrum and is handled by the lag evaluation.
%     .conj      (true)   the stored data are in the engineering sign
%                         (OPR's matched-filter output, the 'deramped'
%                         flag elsewhere): the phase of <a conj(b)>
%                         DEcreases as a is delayed. The data are
%                         conjugated into the physics sign, where it
%                         increases. Set false for data already there.
%     .win       (round(50 m / dz) samples) window length
%     .step      (win/2)  window step
%     .max_lag_ns (8)     lag search half-range
%     .lag_step_ns (0.05) lag grid spacing
%     .t0        (1)      sample index where the delay is known to be zero
%                         (the surface or the firn reference window):
%                         tracking starts from the window containing it.
%     .min_gamma (0.1)    windows below this normalised correlation are
%                         not tracked through (NaN) and do not steer the
%                         fringe integer.
%
% Output
%   out.idx      [Nw x 1] window-centre sample index
%   out.dtau     [Nw x 1] tracked travel-time difference [ns] (NaN where
%                         the correlation is too weak)
%   out.fine     [Nw x 1] the carrier-phase delay in (-fringe/2, fringe/2]
%                         before any fringe integer is chosen: wrapped
%                         differences of it between neighbouring windows
%                         are the local slope, free of any tracking
%   out.dtau_env [Nw x 1] envelope-only estimate [ns], the coarse check
%   out.gamma    [Nw x 1] normalised correlation magnitude at the peak
%   out.sigma    [Nw x 1] 1-sigma delay uncertainty [ns] from the
%                         coherence and the number of independent looks
%   out.n_looks  [Nw x 1] independent looks per window (range cells x
%                         independent traces)
%   out.f_rf     [Nw x 1] the radio frequency [Hz] the pair's ZERO-LAG
%                         correlation turns at: fc minus the centroid of
%                         the cross-spectrum magnitude. The carrier phase
%                         read at the envelope peak counts the delay in
%                         fringes of fc exactly, but anything read at zero
%                         lag - the cross-polarised power 45 deg off the
%                         axes, so the half-wave nodes - counts it in
%                         fringes of the band's own centre, which on the
%                         accumulation radar sits ~13 MHz below the LO.
%   out.power_a, out.power_b [Nw x 1] mean power per sample of each signal
%                         in the window (their ratio is the reflection
%                         anisotropy between the two polarisations)
%   out.env_disagree [Nw x 1] logical: envelope peak > half a fringe from
%                         the tracked delay
%   out.lags     [1 x Nl] the lag grid [ns]
%   out.corr     [Nw x Nl] normalised envelope correlation, the poster's
%                         image
%   out.fringe_ns          one carrier fringe, 1e9 / fc
if nargin < 4, opts = struct(); end
if ~isequal(size(a), size(b))
  error('ptt:xcorrTraveltime:shape', 'a and b must have the same size');
end
fc = H_opt(opts, 'fc', 750e6);
do_conj = H_opt(opts, 'conj', true);
[Nt, Nx] = size(a);
win = H_opt(opts, 'win', max(8, round(50 / (dt * 0.299792458 / (2 * sqrt(3.171))))));
win = min(win, Nt);
step = H_opt(opts, 'step', max(1, round(win / 2)));
max_lag = H_opt(opts, 'max_lag_ns', 8);
lag_step = H_opt(opts, 'lag_step_ns', 0.05);
t0 = H_opt(opts, 't0', 1);
min_gamma = H_opt(opts, 'min_gamma', 0.1);
fringe = 1e9 / fc;

a = double(a); b = double(b);
if do_conj
  % back to the physics sign: a delay advances the phase of <a conj(b)>
  a = conj(a); b = conj(b);
end

lags = -max_lag:lag_step:max_lag;
Nl = numel(lags);
starts = 1:step:max(1, Nt - win + 1);
Nw = numel(starts);
idx = starts(:) + floor(win / 2);

% the window's baseband frequency axis, for lag evaluation in the
% frequency domain. In the physics sign a delayed a is a(t) = b(t - tau)
% exp(+2 pi i fc tau), so A = B exp(-2 pi i f tau) exp(2 pi i fc tau) and
%   C(tau') = sum_f A conj(B) exp(+2 pi i f tau')
% peaks in magnitude at tau' = tau (the envelope), where its phase is the
% carrier term 2 pi fc tau that the mixing left behind (the fine delay,
% modulo a fringe). Evaluating the phase at the envelope peak rather than
% at tau itself leaves a residual ramp of the baseband centroid (-13 MHz)
% times the envelope error: under 0.01 ns for an error of half a fringe.
Nf = 2 ^ nextpow2(2 * win);
f = ((0:Nf-1) - floor(Nf / 2)) / (Nf * dt);      % GHz = cycles per ns
f = ifftshift(f).';
E = exp(2i * pi * f * lags);                      % [Nf x Nl]
hann = 0.5 - 0.5 * cos(2 * pi * (0:win-1).' / (win - 1));
taper_eff = sum(hann) ^ 2 / (win * sum(hann .^ 2));

dtau_env = nan(Nw, 1); gamma = nan(Nw, 1); phase_pk = nan(Nw, 1);
n_looks = nan(Nw, 1); corr = nan(Nw, Nl);
f_x = nan(Nw, 1); power_a = nan(Nw, 1); power_b = nan(Nw, 1);
p_norm = Nf * Nx * sum(hann .^ 2);                % Parseval: sum |A|^2 over bins and traces -> mean power per sample
for w = 1:Nw
  r = starts(w):starts(w) + win - 1;
  A = fft(hann .* a(r, :), Nf, 1);
  B = fft(hann .* b(r, :), Nf, 1);
  X = sum(A .* conj(B), 2);                       % cross-spectrum, traces pooled
  Paa = sum(abs(A(:)) .^ 2); Pbb = sum(abs(B(:)) .^ 2);
  if ~(Paa > 0 && Pbb > 0), continue; end
  power_a(w) = Paa / p_norm; power_b(w) = Pbb / p_norm;
  % the baseband frequency the correlation's phase turns at with lag: the
  % centroid of |X|. Each baseband f carries the radio frequency fc - f in
  % this sign, so the zero-lag correlation winds with delay at fc - f_x
  f_x(w) = sum(f .* abs(X)) / sum(abs(X));
  C = (X.' * E) / sqrt(Paa * Pbb);                % [1 x Nl], normalised
  corr(w, :) = abs(C);
  [g, k] = max(abs(C));
  gamma(w) = g;
  % parabolic refinement of the envelope peak
  if k > 1 && k < Nl
    y0 = abs(C(k-1)); y1 = abs(C(k)); y2 = abs(C(k+1));
    den = y0 - 2 * y1 + y2;
    sh = 0; if den < 0, sh = 0.5 * (y0 - y2) / den; end
    dtau_env(w) = lags(k) + sh * lag_step;
  else
    dtau_env(w) = lags(k);
  end
  % carrier phase at the envelope peak: 2 pi fc tau modulo one fringe
  phase_pk(w) = angle(C(k));
  % independent looks: range cells (effective window length x effective
  % bandwidth, the Hann taper shortening the window by sum(h)^2 / (N sum(h^2)))
  % times independent traces (along-track correlation of a)
  Pf = mean(abs(A) .^ 2, 2);
  bw = sum(Pf) ^ 2 / sum(Pf .^ 2) / (Nf * dt);    % effective bandwidth, GHz
  nx_eff = Nx / max(1, H_tau_x(a(r, :)));
  n_looks(w) = max(1, win * taper_eff * dt * bw) * nx_eff;
end

% fine delay: the carrier phase in fringes, one branch of many; the fringe
% integer is tracked by continuity from the zero-delay reference and the
% envelope lag is kept as the independent check
fine = phase_pk / (2 * pi) * fringe;
dtau = nan(Nw, 1);
[~, w0] = min(abs(idx - t0));
prev = 0;
order = [w0:Nw, w0-1:-1:1];
for w = order
  if w == w0 - 1, prev = 0; end
  if ~(gamma(w) >= min_gamma), continue; end
  m = round((prev - fine(w)) / fringe);
  dtau(w) = fine(w) + m * fringe;
  prev = dtau(w);
end
env_disagree = abs(dtau - dtau_env) > fringe / 2;
% delay uncertainty from the phase: Cramer-Rao for a coherence gamma over
% n independent looks, sigma_phi = sqrt(1 - g^2) / (g sqrt(2 n))
sigma = sqrt(max(0, 1 - gamma .^ 2)) ./ max(gamma, 1e-6) ./ sqrt(2 * n_looks) ...
  / (2 * pi) * fringe;

out = struct('idx', idx, 'dtau', dtau, 'fine', fine, 'dtau_env', dtau_env, 'gamma', gamma, ...
  'sigma', sigma, 'n_looks', n_looks, 'env_disagree', env_disagree, ...
  'f_rf', fc - f_x * 1e9, 'power_a', power_a, 'power_b', power_b, ...
  'lags', lags, 'corr', corr, 'fringe_ns', fringe, 'win', win, 'step', step);
end

function v = H_opt(o, name, default)
if isfield(o, name) && ~isempty(o.(name)), v = o.(name); else, v = default; end
end

function tau = H_tau_x(S)
% along-track correlation length in traces, as the integrated squared
% correlation 1 + 2 sum_L |rho(L)|^2 - the factor by which correlated
% speckle reduces the number of independent traces in a variance (for
% Gaussian speckle the variance of a cross-moment scales with |rho|^2)
Nx = size(S, 2);
if Nx < 4, tau = 1; return; end
m = sum(abs(S) .^ 2, 1);
r = zeros(1, min(Nx - 1, 50));
for L = 1:numel(r)
  r(L) = abs(sum(sum(S(:, 1:end-L) .* conj(S(:, 1+L:end))))) / sqrt(sum(m(1:end-L)) * sum(m(1+L:end)));
end
% stop the sum where the correlation has fallen into its noise floor
k = find(r < 0.1, 1);
if ~isempty(k), r = r(1:k); end
tau = 1 + 2 * sum(r .^ 2);
end
