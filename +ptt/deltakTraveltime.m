function [dtau, info] = deltakTraveltime(slc, map, opts)
%DELTAKTRAVELTIME Absolute polarimetric dtau by split-spectrum delta-k.
%   [dtau, info] = DELTAKTRAVELTIME(slc, map, opts) estimates the two-way
%   traveltime difference dtau(fast-time, along-track) [s] between two
%   polarization channels WITHOUT phase unwrapping, from the linearity of
%   the interferometric phase in RF frequency (the delay is non-
%   dispersive): phi(f) = s*2*pi*f*dtau. Sub-band interferograms are
%   cross-multiplied PER PIXEL (their difference phase varies as
%   2*pi*df*dtau - slow everywhere, even where full-carrier fringes are
%   dense) and only then multilooked, in a coarse-to-fine ladder:
%     stage A: adjacent-pair cross products over n_sub narrow sub-bands
%              (df ~ B/n_sub; unambiguous |dtau| < 1/(2*df))
%     stage Q: quarter-band pairs, integers resolved by smoothed A
%     stage B: half-band pair, integers resolved by Q -> absolute dtau
%   The result is referenced to zero over a shallow band below the
%   surface (removes channel timing/phase bias) and returned interpolated
%   to the full (Nt x Nx) grid for ptt.blockAverage.
%
%   slc fields:
%     ref   complex SLC of the reference channel (Nt x Nx)
%     sec   complex SLC of the secondary channel, UNREGISTERED: image
%           coregistration shifts the envelope, which removes exactly the
%           group delay this method measures. Narrow sub-bands tolerate
%           the few-bin misregistration (envelope correlation width
%           1/B_sub is much larger than the offsets).
%
%   map: as for ptt.blendTraveltime (Time, Surface, fc, coherence, and
%   optionally row_offset for the orientation regression).
%
%   opts fields used (all optional, under opts.deltak unless noted):
%     n_sub (12)            narrow sub-bands for stage A
%     cell_twtt (100e-9)    analysis-cell fast-time extent [s]
%     cell_ntr (31)         analysis-cell trace count
%     smooth (5)            resolver smoothing (cells, boxcar)
%     band_frac ([.02 .98]) cumulative-power band support
%     orientation (0)       spectral orientation s; 0 = regress against
%                           map.row_offset (row_offset > 0 = sec later)
%     ref_band ([150e-9 600e-9])  referencing band below the surface
%     coherence_threshold (opts.*, 0.5)  as in blendTraveltime
%
%   info: phase_sign (the detected/forced orientation s), coh_mask,
%   dtau_coreg = [] (delta-k needs no fringe blending), phase_is_unwrapped
%   = true, ref_bin, and diagnostics under info.deltak (cell-grid taus,
%   band frequencies, ladder rounding margins).
%
%   The full-carrier refinement (resolving the fc fringe integer from
%   stage B) is deliberately NOT applied: on real data the carrier phase
%   shows a residual inconsistency with the group delay at the ~1/4
%   fringe level (registration application details / channel dispersion),
%   so the synthetic-wavelength stage-B product is the deliverable.

if ~isfield(opts,'deltak') || isempty(opts.deltak)
  opts.deltak = struct();
end
dk = opts.deltak;
if ~isfield(dk,'n_sub') || isempty(dk.n_sub), dk.n_sub = 12; end
if ~isfield(dk,'cell_twtt') || isempty(dk.cell_twtt), dk.cell_twtt = 100e-9; end
if ~isfield(dk,'cell_ntr') || isempty(dk.cell_ntr), dk.cell_ntr = 31; end
if ~isfield(dk,'smooth') || isempty(dk.smooth), dk.smooth = 5; end
if ~isfield(dk,'band_frac') || isempty(dk.band_frac), dk.band_frac = [0.02 0.98]; end
if ~isfield(dk,'orientation') || isempty(dk.orientation), dk.orientation = 0; end
if ~isfield(dk,'ref_band') || isempty(dk.ref_band), dk.ref_band = [150e-9 600e-9]; end
if ~isfield(opts,'coherence_threshold') || isempty(opts.coherence_threshold)
  opts.coherence_threshold = 0.5;
end
if ~isfield(opts,'ref_twtt_offset') || isempty(opts.ref_twtt_offset)
  opts.ref_twtt_offset = 50e-9;
end

[Nt, Nx] = size(slc.ref);
dt = map.Time(2) - map.Time(1);
fs = 1/dt;
fc = map.fc;

%% Range spectra and band support
R = fft(single(slc.ref), [], 1);
S = fft(single(slc.sec), [], 1);
pw = mean(abs(R(:, 1:max(1,floor(Nx/200)):end)).^2, 2);
pw = fftshift(pw);
fbb = ((0:Nt-1).' - floor(Nt/2)) * fs / Nt;
cum = cumsum(pw) / sum(pw);
b0 = find(cum >= dk.band_frac(1), 1);
b1 = find(cum >= dk.band_frac(2), 1);
fprintf('deltak: band support %.0f..%.0f MHz baseband (%.0f MHz)\n', ...
  fbb(b0)/1e6, fbb(b1)/1e6, (fbb(b1)-fbb(b0))/1e6);

tg = max(1, round(dk.cell_twtt/dt));  % cell size in fast-time samples
xg = dk.cell_ntr;
ntc = floor(Nt/tg);
nxc = floor(Nx/xg);
t_cell = map.Time(1) + ((0:ntc-1).' + 0.5)*tg*dt;
x_cell = ((0:nxc-1) + 0.5)*xg;

%% Ladder stages: per-pixel cross products, multilook after differencing
edgesA = round(linspace(b0, b1, dk.n_sub+1));
[accA, dfA] = stage_cross(edgesA);
edgesQ = round(linspace(b0, b1, 5));
[accQ, dfQ] = stage_cross(edgesQ);
edgesB = round(linspace(b0, b1, 3));
[accB, dfB] = stage_cross(edgesB);
fprintf('deltak: df A/Q/B = %.1f/%.1f/%.1f MHz (unambiguous +/- %.0f/%.0f/%.1f ns)\n', ...
  dfA/1e6, dfQ/1e6, dfB/1e6, 1e9/(2*dfA), 1e9/(2*dfQ), 1e9/(2*dfB));

k = ones(dk.smooth)/dk.smooth^2;
tau_A = angle(conv2(real(accA), k, 'same') + 1i*conv2(imag(accA), k, 'same')) ...
  / (2*pi*dfA);
tau_Qw = angle(conv2(real(accQ), k, 'same') + 1i*conv2(imag(accQ), k, 'same')) ...
  / (2*pi*dfQ);
k3 = ones(3)/9;
tau_Bw = angle(conv2(real(accB), k3, 'same') + 1i*conv2(imag(accB), k3, 'same')) ...
  / (2*pi*dfB);

%% Spectral orientation: regression of the group delay against the
% coregistration offsets (row_offset > 0 = sec arrives later)
coh_c = cellavg(map.coherence);
surf_c = cellavg_row(map.Surface);
if dk.orientation ~= 0
  s = sign(dk.orientation);
elseif isfield(map,'row_offset') && ~isempty(map.row_offset)
  coreg_c = cellavg(map.row_offset*dt);
  tau_A_r = band_ref(tau_A);
  coreg_r = band_ref(coreg_c);
  msk = coh_c > 0.35 & isfinite(coreg_r) & abs(coreg_r) > dt/2;
  if nnz(msk) > 100
    s = sign(sum(tau_A_r(msk).*coreg_r(msk)));
    if s == 0, s = -1; end
  else
    warning('ptt:deltakTraveltime:sign', ...
      'Too few coregistration cells for orientation; using -1 (matched filter).');
    s = -1;
  end
else
  warning('ptt:deltakTraveltime:sign', ...
    'No row_offset for orientation regression; using -1 (matched filter).');
  s = -1;
end
tau_A = s*tau_A; tau_Qw = s*tau_Qw; tau_Bw = s*tau_Bw;

%% Resolve the ladder integers per cell
nQ = round((tau_A - tau_Qw)*dfQ);
tau_Q = tau_Qw + nQ/dfQ;
nB = round((tau_Q - tau_Bw)*dfB);
tau_B = tau_Bw + nB/dfB;
mQ = abs((tau_A - tau_Qw)*dfQ - nQ);
mB = abs((tau_Q - tau_Bw)*dfB - nB);
fprintf('deltak: rounding margins Q %.2f/%.2f, B %.2f/%.2f (median/90th; <<0.5)\n', ...
  nmed(mQ(:)), prctile_ish(mQ(:),90), nmed(mB(:)), prctile_ish(mB(:),90));

%% Surface referencing (robust shallow-band median per cell column)
tau_ref = band_ref(tau_B);

%% Interpolate to the full grid for ptt.blockAverage
ti = min(max(interp1(t_cell, (1:ntc).', map.Time(:), 'linear', 'extrap'), 1), ntc);
xi = min(max(interp1(x_cell, (1:nxc).', (1:Nx).', 'linear', 'extrap'), 1), nxc);
dtau = interp2(tau_ref, xi.', ti, 'linear');

surf_valid = isfinite(map.Surface(:).');
coh_mask = map.coherence >= opts.coherence_threshold;
coh_mask(:,~surf_valid) = false;
ref_bin = round((map.Surface(:).' + opts.ref_twtt_offset - map.Time(1))/dt) + 1;
ref_bin(~surf_valid) = 1;
ref_bin = min(max(ref_bin,1),Nt);

info.phase_sign = s;
info.coh_mask = coh_mask;
info.dtau_coreg = [];   % no fringe blending: dtau is already absolute
info.phase_is_unwrapped = true;
info.ref_bin = ref_bin;
info.deltak = struct('tau_cell', tau_ref, 't_cell', t_cell, ...
  'x_cell', x_cell, 'df', [dfA dfQ dfB], ...
  'margin_Q', [nmed(mQ(:)) prctile_ish(mQ(:),90)], ...
  'margin_B', [nmed(mB(:)) prctile_ish(mB(:),90)]);

%% ---- nested helpers ------------------------------------------------
  function [acc, df] = stage_cross(edges)
    % Per-pixel cross products of adjacent band-pass interferograms,
    % accumulated, then multilooked onto the cell grid.
    acc_full = [];
    prev = [];
    fcs = zeros(1, numel(edges)-1);
    for m = 1:numel(edges)-1
      sl = edges(m):edges(m+1)-1;
      un = mod(sl - 1 - floor(Nt/2), Nt) + 1;  % shifted -> unshifted bins
      Rb = complex(zeros(Nt, Nx, 'single'));
      Sb = complex(zeros(Nt, Nx, 'single'));
      Rb(un,:) = R(un,:);
      Sb(un,:) = S(un,:);
      Ib = ifft(Sb, [], 1) .* conj(ifft(Rb, [], 1));
      fcs(m) = fc + sum(fbb(sl).*pw(sl))/sum(pw(sl));
      if ~isempty(prev)
        x = Ib .* conj(prev);
        if isempty(acc_full), acc_full = x; else, acc_full = acc_full + x; end
      end
      prev = Ib;
    end
    acc = cellavg(acc_full);
    df = mean(diff(fcs));
  end

  function C = cellavg(A)
    % Mean over (tg x xg) cells; A is (Nt x Nx) -> (ntc x nxc)
    B = reshape(A(1:ntc*tg, 1:nxc*xg), tg, ntc, xg, nxc);
    C = squeeze(mean(mean(B, 1), 3));
    C = reshape(C, ntc, nxc);   % squeeze safety for ntc==1 or nxc==1
  end

  function m = nmed(v)
    % NaN-tolerant median (Octave-safe: no 'omitnan' dependence)
    v = v(isfinite(v));
    if isempty(v), m = NaN; else, m = median(v); end
  end

  function r = cellavg_row(v)
    % Cell means of a per-trace row vector (1 x Nx) -> (1 x nxc)
    r = mean(reshape(v(1:nxc*xg), xg, nxc), 1);
  end

  function A = band_ref(A)
    % Subtract a robust shallow-band median per cell column
    for i = 1:nxc
      band = t_cell > surf_c(i) + dk.ref_band(1) & ...
        t_cell < surf_c(i) + dk.ref_band(2);
      if any(band)
        A(:,i) = A(:,i) - nmed(A(band,i));
      end
    end
  end

  function p = prctile_ish(v, q)
    % Percentile without the stats toolbox (Octave-safe)
    v = sort(v(isfinite(v)));
    if isempty(v), p = NaN; return; end
    p = v(max(1, min(numel(v), round(q/100*numel(v)))));
  end

end
