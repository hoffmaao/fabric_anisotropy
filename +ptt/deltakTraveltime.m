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
%   opts fields used (all optional, under opts.deltak unless noted, and
%   defaulted by ptt.deltakDefaults so ptt.imgCombSeam can widen the seam
%   band by the same analysis-cell reach):
%     n_sub (12)            narrow sub-bands for stage A
%     cell_twtt (100e-9)    analysis-cell fast-time extent [s]
%     cell_ntr (31)         analysis-cell trace count
%     smooth (5)            resolver smoothing (cells, boxcar)
%     band_frac ([.02 .98]) cumulative-power band support
%     orientation (0)       spectral orientation s; 0 = regress against
%                           map.row_offset (row_offset > 0 = sec later)
%     phase_sign (top-level opts.phase_sign, 0) forced orientation
%                           override, as in ptt.blendTraveltime; nonzero
%                           takes precedence over orientation and the
%                           regression
%     ref_band ([150e-9 600e-9])  referencing band below the surface
%     ref_min_cells (2)     analysis cells the referencing band must keep
%                           after the waveform-combine seam cells are
%                           excluded; below it the column is unreferenced
%     ref_extend_max (1e-6) how much deeper than ref_band(2) the band may
%                           step to reach ref_min_cells usable cells when
%                           the seam ate into it; past the cap the column
%                           is unreferenced rather than referenced off
%                           something that is no longer a shallow band
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

dk = ptt.deltakDefaults(opts);
% coherence_threshold and ref_twtt_offset are defaulted by
% ptt.surfaceReference, which owns the referencing convention

[Nt, Nx] = size(slc.ref);
dt = map.Time(2) - map.Time(1);
fs = 1/dt;
fc = map.fc;

%% Range spectra and band support
% Non-finite samples must be zeroed first: the FFT would otherwise spread a
% single NaN over that whole trace, and the 5x5 cell smoothing below over
% its neighbours, voiding an entire block average downstream.
R = single(slc.ref);
S = single(slc.sec);
bad = ~isfinite(R) | ~isfinite(S);
if any(bad(:))
  warning('ptt:deltakTraveltime:nonfinite', ...
    'Zeroing %d non-finite SLC samples before the range FFT.', nnz(bad));
  R(bad) = 0;
  S(bad) = 0;
end
bad = [];
R = fft(R, [], 1);
S = fft(S, [], 1);
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

% The ladder integers come from a tau_A smoothed over the analysis cells,
% so a cell whose smoothing window merely touched a waveform-combine seam
% is already resolved to the wrong 1/dfQ step: ask for the widened band
% (ptt.imgCombSeam takes the reach from the same ptt.deltakDefaults). The
% same cells are kept out of the referencing median in band_ref, so the
% estimator and ptt.surfaceReference agree on which cells are usable.
opts.seam_mask_deltak = true;
seam_cell = seam_cells();
n_ref_extended = 0;

%% Ladder stages: per-pixel cross products, multilook after differencing
edgesA = round(linspace(b0, b1, dk.n_sub+1));
[accA, dfA] = stage_cross(edgesA);
edgesQ = round(linspace(b0, b1, 5));
[accQ, dfQ] = stage_cross(edgesQ);
edgesB = round(linspace(b0, b1, 3));
[accB, dfB] = stage_cross(edgesB);
R = []; S = [];   % full-frame spectra are done with; release before smoothing
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
if isfield(opts,'phase_sign') && ~isempty(opts.phase_sign) && opts.phase_sign ~= 0
  s = sign(opts.phase_sign);
elseif dk.orientation ~= 0
  s = sign(dk.orientation);
elseif isfield(map,'row_offset') && ~isempty(map.row_offset)
  coreg_c = cellavg(map.row_offset*dt);
  tau_A_r = band_ref(tau_A);
  coreg_r = band_ref(coreg_c);
  msk = coh_c > 0.35 & isfinite(coreg_r) & isfinite(tau_A_r) ...
    & abs(coreg_r) > dt/2;
  if nnz(msk) > 100
    s = sign(sum(tau_A_r(msk).*coreg_r(msk)));
    if ~isfinite(s) || s == 0
      warning('ptt:deltakTraveltime:sign', ...
        'Orientation regression is degenerate; using -1 (matched filter).');
      s = -1;
    elseif s ~= -1
      warning('ptt:deltakTraveltime:sign', ...
        'Orientation regression detected s = %+d, contrary to the matched-filter convention (-1); set opts.phase_sign to force.', s);
    end
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
if n_ref_extended > 0
  warning('ptt:deltakTraveltime:refExtended', ...
    '%d of %d analysis-cell columns kept fewer than %g reference cells inside %.0f ns of the surface, so the band was stepped up to %.0f ns deeper to reach them.', ...
    n_ref_extended, nxc, dk.ref_min_cells, dk.ref_band(2)*1e9, ...
    dk.ref_extend_max*1e9);
end
n_unref = nnz(all(~isfinite(tau_ref), 1));
if n_unref > 0
  warning('ptt:deltakTraveltime:unreferenced', ...
    '%d of %d analysis-cell columns keep fewer than %d usable reference cells (no surface pick, band outside the record, or the waveform-combine seam) and are excluded from dtau.', ...
    n_unref, nxc, dk.ref_min_cells);
end

%% Interpolate to the full grid for ptt.blockAverage
ti = min(max(interp1(t_cell, (1:ntc).', map.Time(:), 'linear', 'extrap'), 1), ntc);
xi = min(max(interp1(x_cell, (1:nxc).', (1:Nx).', 'linear', 'extrap'), 1), nxc);
dtau = interp2(tau_ref, xi.', ti, 'linear');

[coh_mask, ref_bin] = ptt.surfaceReference(map, opts);

% Cells left unreferenced (no surface pick anywhere in the cell) stay
% non-finite: drop them from the mask and zero them so that a single bad
% cell cannot NaN a whole block average in ptt.blockAverage.
bad_dtau = ~isfinite(dtau);
coh_mask(bad_dtau) = false;
dtau(bad_dtau) = 0;

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
    % One pair of band buffers for the whole ladder stage: only the bins of
    % the band in hand are ever non-zero, so they are cleared again after
    % each band rather than reallocated (two full-frame arrays per band).
    Rb = complex(zeros(Nt, Nx, 'single'));
    Sb = complex(zeros(Nt, Nx, 'single'));
    for m = 1:numel(edges)-1
      sl = edges(m):edges(m+1)-1;
      un = mod(sl - 1 - floor(Nt/2), Nt) + 1;  % shifted -> unshifted bins
      Rb(un,:) = R(un,:);
      Sb(un,:) = S(un,:);
      Ib = ifft(Sb, [], 1) .* conj(ifft(Rb, [], 1));
      Rb(un,:) = complex(0,0);
      Sb(un,:) = complex(0,0);
      fcs(m) = fc + sum(fbb(sl).*pw(sl))/sum(pw(sl));
      if ~isempty(prev)
        x = Ib .* conj(prev);
        if isempty(acc_full), acc_full = x; else, acc_full = acc_full + x; end
      end
      prev = Ib;
    end
    acc = cellavg(acc_full);
    df = mean(diff(fcs));
    % Nested-function variables live in the parent workspace, so the
    % full-frame temporaries would otherwise stay resident across stages
    Rb = []; Sb = []; Ib = []; prev = []; x = []; acc_full = [];
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
    % Cell medians of a per-trace row vector (1 x Nx) -> (1 x nxc). NaN
    % tolerant: one bad pick must not void the whole cell's referencing.
    Bv = reshape(v(1:nxc*xg), xg, nxc);
    r = nan(1, nxc);
    for ic = 1:nxc
      r(ic) = nmed(Bv(:,ic));
    end
  end

  function C = seam_cells()
    % Analysis cells the widened waveform-combine seam band covers, as
    % (ntc x 1) when the band is column-independent and (ntc x nxc) when a
    % boundary tracks the surface. The band itself is never recomputed
    % here: it comes from ptt.imgCombSeam, the same mask
    % ptt.surfaceReference applies to the coherence. A cell is judged by
    % the sample at its centre, t_cell, because the widened band already
    % carries the neighbourhood reach - testing every sample in the cell
    % would count that half-cell twice and swallow the cell next door.
    if isfield(opts,'seam_mask_en') && ~isempty(opts.seam_mask_en) ...
        && ~opts.seam_mask_en
      C = false(ntc, 1);
      return;
    end
    sm = ptt.imgCombSeam(map, opts);
    if ~any(sm(:))
      C = false(ntc, 1);
      return;
    end
    ci = (0:ntc-1).'*tg + floor(tg/2) + 1;   % bin at each t_cell
    if size(sm,2) == 1
      C = sm(ci);
    else
      % ptt.imgCombSeam resolves the surface-tracking case per TRACE;
      % reduce onto the cell columns band_ref indexes, the way cellavg
      % and cellavg_row do, so a cell column is judged by its own traces
      Bs = reshape(sm(ci, 1:nxc*xg), ntc, xg, nxc);
      C = reshape(any(Bs, 2), ntc, nxc);
      Bs = [];
    end
    sm = [];
  end

  function A = band_ref(A)
    % Subtract a robust shallow-band median per cell column, over the cells
    % the seam mask leaves usable. Where the seam has eaten into the
    % declared band the far edge steps DEEPER - never shallower, the near
    % edge is pinned by the surface - by at most dk.ref_extend_max, taking
    % the shallowest dk.ref_min_cells usable cells it finds. That keeps a
    % column referenced instead of lost when the record starts too late for
    % the declared window to hold enough cells on its own; the cap is there
    % because this is supposed to be a shallow-band zero level, and a
    % reference taken that much deeper already carries whatever dtau the
    % column has accumulated by then.
    %
    % A column that still keeps fewer than dk.ref_min_cells (no finite
    % Surface pick in the cell, the band falling outside the record, or the
    % seam swallowing it beyond the cap) is NaN rather than unreferenced:
    % the channel timing bias it would otherwise keep is not a delay, and
    % interp2 would smear it into the valid cells on either side. That is
    % per column - it costs those traces, never the frame - and the count
    % is reported by the ptt:deltakTraveltime:unreferenced warning.
    % Referencing off a single surviving cell is not robust enough to be
    % worth preferring to either.
    nmin = max(1, round(dk.ref_min_cells));
    n_ref_extended = 0;
    for i = 1:nxc
      if size(seam_cell,2) == 1
        usable = ~seam_cell;
      else
        usable = ~seam_cell(:,i);
      end
      usable = usable & t_cell > surf_c(i) + dk.ref_band(1);
      band = usable & t_cell < surf_c(i) + dk.ref_band(2);
      if nnz(band) < nmin
        deeper = find(usable & ...
          t_cell < surf_c(i) + dk.ref_band(2) + dk.ref_extend_max);
        if numel(deeper) >= nmin
          band = false(ntc,1);
          band(deeper(1:nmin)) = true;
          n_ref_extended = n_ref_extended + 1;
        end
      end
      if nnz(band) >= nmin
        A(:,i) = A(:,i) - nmed(A(band,i));
      else
        A(:,i) = NaN;
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
