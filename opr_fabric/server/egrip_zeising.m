%EGRIP_ZEISING Zeising et al. (2023) phase co-registration fabric estimator.
%
% Reference implementation of the METHOD WE ARE BEING COMPARED AGAINST, run
% on our own EastGRIP lines so the comparison is like-for-like. Zeising, O.,
% Gerber, T. A., Eisen, O., Ershadi, M. R., Stoll, N., Weikusat, I. and
% Humbert, A. (2023), "Improved estimation of the bulk ice crystal fabric
% asymmetry from polarimetric phase co-registration", The Cryosphere 17,
% 1097-1105, doi:10.5194/tc-17-1097-2023.
%
% Their estimator, eq. (2)-(4) of the paper, is
%
%   c_hhvv(z,l) = sum_j s_hh(j) conj(s_vv(j+l))
%                 / sqrt(sum |s_hh|^2 * sum |s_vv|^2)
%   dt(z)       = l(z)/(B*p) + arg(c_hhvv)/(2*pi*fc)
%   dlam(z)     = (d dt/dz) * c * n / deps'        with deps' = 0.034
%
% i.e. a coarse bin lag that resolves which whole cycle the pair sits in,
% refined by the coherence phase inside that cycle. That is structurally the
% same combination ptt.blendTraveltime already makes at
%   fringe = round((dtau_coreg - dtau_phase)*fc)
% so this script is NOT a different physical model - it is the same delay,
% estimated pointwise and differentiated locally, where our own path instead
% inverts the whole dtau(z) profile under a regulariser and can additionally
% use the returned power. Keeping the physics identical is what makes the
% remaining difference attributable to the inversion.
%
% WHY THIS RUNS ON THE SERVER. The lagged product s_hh(j) conj(s_vv(j+l))
% cannot be formed from the staged extracts: those carry the ZERO-LAG
% product interferogram_mlook = sum HH conj(VV) and the two power sums, and
% no lagged product is recoverable from them. The estimator needs the
% separate complex channels at FULL range resolution, so it reads
% CSARP_qlook_HH/VV directly and returns only the small per-line profiles.
%
% AZIMUTH AGGREGATION. Their sounding is stationary; ours moves. The sum
% over j is extended over the near-borehole traces as well, which is valid
% for the same reason the sum over j is: at the correct lag the layer phase
% common to HH and VV cancels inside the product, so each trace's segment
% sum carries the same sub-bin residual phase and they add coherently. This
% is NOT true of a coherent stack of the channels themselves - two traces a
% wavelength apart (0.224 m in ice at 750 MHz) carry independent layer
% phase, so the channels are never stacked before the product is formed.
%
% WHERE THE SUB-BIN LAG COMES FROM. Their padding factor p=8 buys a 0.625 ns
% bin against a 3.33 ns cycle at 300 MHz, so their coarse lag alone says
% which cycle the pair is in. Our raw bin is 1.667 ns against a 1.333 ns
% cycle at 750 MHz - the opposite regime, where an unpadded lag cannot. The
% margin is restored by interpolating the CORRELATION rather than the data:
% padding the profiles x8 before correlating would cost ~3 GB per line and
% is mathematically the same thing, since correlation and band-limited
% interpolation commute. The correlation is 43 lags long, so interpft on it
% is free.
%
% Positions, the stopped-trace cull and the surface pick are all identical
% to negis_interferogram.m; see that file for why each is needed.
season = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
out_dir = '/kucresis/scratch/hoffmana_sta/fabric/stages';

% The nine borehole-proximal lines, and they must stay in step with FRAMES
% in scripts/figures/egrip_azimuthal.py: the azimuthal solve that consumes
% this output fits a two-parameter cos2a/sin2a model over exactly these
% lines, and eight of the nine sit between 128 and 145 deg, so the single
% line near 33 deg carries all of the orthogonal leverage.
frames = { '20240628_01', 1; '20240626_03', 5; '20240626_01', 1; ...
           '20240619_01', 1; '20240621_01', 1; '20240620_01', 1; ...
           '20240621_01', 10; '20240619_01', 5; '20240620_02', 3 };

EG_LAT = 75.6294; EG_LON = -35.9672;   % camp survey centroid (20240618_01)
RADIUS_KM = 2.0;        % must match RADIUS_KM in egrip_azimuthal.py
MIN_STEP_M = 1.0;       % below this along-track step the vehicle was stopped
MAX_TRACES = 1500;      % decimate above this; SNR is already saturated

% --- Zeising et al. (2023) parameters, section 3
% Segment length is the one parameter that CANNOT be carried over as a
% length. Theirs is 12 m at fc = 300 MHz; the quantity that matters is how
% much the sub-bin phase drifts ACROSS a segment, and that scales with fc.
% At the core's dlam ~ 0.5 a 12 m segment drifts 0.29 cycles (103 deg) at
% our 750 MHz against 0.115 cycles (41 deg) at theirs, which both lowers
% |c| and biases its angle. 12*300/750 = 4.8 m reproduces their phase
% budget. All three are run and reported so the choice is visible rather
% than asserted - and so "their method scored X" cannot be an artefact of
% handing it a segment length their radar never implied.
SEG_LIST  = [12.0, 6.0, 4.8];
COH_MIN   = 0.65;       % their correlation threshold
SMOOTH_M  = 100;        % moving average applied to dt
DERIV_M   = 200;        % moving window for the depth derivative of dt
PAD       = 8;          % their padding factor p, here on the correlation
DEPS      = 0.034;      % dielectric anisotropy of the ice crystal
% Lag half-width. The column accumulates dt = 6.39e-11 * dlam * z, so at
% the core's dlam ~ 0.5 over 1400 m that is 44.7 ns = 3.8 m of ice. A +-3 m
% window (the first cut here) saturates below ~1100 m and silently pins the
% peak to the window edge, which reads as a collapse of dlam at depth.
MAXLAG_M  = 8.0;
GAP_M     = 24;         % re-seed the cycle track across a gap wider than this
Z_TOP     = 12;         % shallowest segment top
Z_BOT     = 1400;       % their deepest core comparison depth

FC = 750e6;             % our centre frequency (theirs is 300 MHz)
C0 = 299792458;
EPS_ICE = 3.15;
N_ICE = sqrt(EPS_ICE);
C_ICE = C0 / N_ICE;

% dlam from the depth gradient of dt. Two-way delay difference accumulates
% at d(dt)/dz = 2*dn/c with dn = deps'*dlam/(2n), so
%   dlam = (d dt/dz) * c * n / deps'.
% Cross-check against the independent calibration the Python figures use
% (K_CONV = 4.021 cycles/us per unit dlam at 750 MHz, measured by calling
% ptt.twttDifference with a uniform dlam in fringe_check.m): that constant
% implies d(dt)/dt_twtt = K_CONV/FC, and dz/dt_twtt = C_ICE/2, so
% d(dt)/dz = 2*K_CONV/(FC*C_ICE) = 6.35e-11 per unit dlam, against
% 2*deps'/(2n*c) = 6.39e-11 here. They agree to 0.7%, so the two methods
% are converted to eigenvalues on the same physics and the comparison is of
% estimators, not of constants.
GRAD_PER_DLAM = 2 * (DEPS / (2 * N_ICE)) / C0;

Rearth = 6371e3;
L = struct([]);
nkeep = 0;

for fi = 1:size(frames, 1)
  day_seg = frames{fi, 1};
  frm = frames{fi, 2};
  name = sprintf('Data_%s_%03d.mat', day_seg, frm);
  tag = sprintf('%s_%03d', day_seg, frm);
  hh_fn = fullfile(season, 'CSARP_qlook_HH', day_seg, name);
  vv_fn = fullfile(season, 'CSARP_qlook_VV', day_seg, name);
  fprintf('\n=== %s ===\n', tag);
  try
    if exist(hh_fn,'file') ~= 2 || exist(vv_fn,'file') ~= 2
      error('egrip_zeising:missing', 'missing HH or VV for %s', name);
    end

    H = load(hh_fn, 'Data', 'Time', 'GPS_time');
    V = load(vv_fn, 'Data', 'Time', 'GPS_time');
    if ~isequal(H.Time, V.Time)
      error('egrip_zeising:timeAxis', 'HH/VV range axes differ');
    end
    if isreal(H.Data) || isreal(V.Data)
      error('egrip_zeising:realData', ...
        'Data is real - this frame was written with incoherent decimation');
    end
    [Nt, Nx] = size(H.Data);
    dt_s = H.Time(2) - H.Time(1);

    % --- reposition from the day GPS file (records time sync failed)
    g = load(fullfile(gps_dir, sprintf('gps_%s.mat', day_seg(1:8))), ...
      'sync_gps_time', 'sync_lat', 'sync_lon');
    [tu, iu] = unique(g.sync_gps_time(:));
    gt = H.GPS_time(:);
    Latitude = interp1(tu, g.sync_lat(iu), gt, 'linear', NaN).';
    Longitude = interp1(tu, g.sync_lon(iu), gt, 'linear', NaN).';

    % --- cull the stops, before anything is averaged
    pp = deg2rad(Latitude); dl = diff(deg2rad(Longitude)); dp = diff(pp);
    aa = sin(dp/2).^2 + cos(pp(1:end-1)).*cos(pp(2:end)).*sin(dl/2).^2;
    step = [Inf, 2*Rearth*asin(sqrt(aa))];
    located = isfinite(Latitude) & isfinite(Longitude);
    stopped = located & [false, located(1:end-1)] & step < MIN_STEP_M;
    keep = located & ~stopped;

    % --- near-borehole aperture
    pb = deg2rad(EG_LAT);
    pl = deg2rad(Latitude); dlo = deg2rad(Longitude - EG_LON);
    ah = sin((pl-pb)/2).^2 + cos(pb).*cos(pl).*sin(dlo/2).^2;
    dist_km = 2*Rearth*asin(sqrt(ah)) / 1e3;
    near = find(keep & (dist_km <= RADIUS_KM));
    fprintf('  %d samples x %d traces; %d near (<= %.1f km), %d culled\n', ...
      Nt, Nx, numel(near), RADIUS_KM, nnz(~keep));
    if numel(near) < 8
      error('egrip_zeising:tooFewNear', ...
        'only %d traces within %.1f km of the borehole', numel(near), ...
        RADIUS_KM);
    end
    if numel(near) > MAX_TRACES
      near = near(round(linspace(1, numel(near), MAX_TRACES)));
    end

    hh = single(H.Data(:, near)); vv = single(V.Data(:, near));
    lat_n = Latitude(near); lon_n = Longitude(near);
    Tvec = H.Time(:);          % keep the real range axis; Time(1) ~= 0
    clear H V;

    % --- surface pick: leading edge, first sample above 5% of trace peak.
    % Identical to negis_interferogram.m and run_negis_fabric.m so this
    % method and ours share a depth zero; a different pick would offset the
    % two dlam(z) curves by the pulse rise and manufacture a disagreement.
    pw = abs(hh).^2;
    thr = 0.05 * max(pw, [], 1);
    surf_s = nan(1, size(hh,2));
    for k = 1:size(hh,2)
      j = find(pw(:,k) > thr(k), 1);
      if ~isempty(j), surf_s(k) = Tvec(j); end
    end
    t_surf = median(surf_s, 'omitnan');
    clear pw thr;

    % Depth of every raw sample, relative to the picked surface
    z_raw = (Tvec - t_surf) * C_ICE / 2;

    maxlag = ceil(MAXLAG_M / (C_ICE/2 * dt_s));
    lags = (-maxlag:maxlag).';
    nlag = numel(lags);

    % --- lagged products, accumulated over traces once per lag.
    % Forming these as full-profile vectors and taking segment sums off a
    % cumulative sum is what makes this tractable: the direct
    % (segment x lag x trace x sample) loop is ~1e13 operations per line.
    Ehh = sum(abs(hh).^2, 2);
    Evv = sum(abs(vv).^2, 2);
    Pl = complex(zeros(Nt, nlag, 'single'));
    for li = 1:nlag
      l = lags(li);
      if l >= 0
        j0 = 1; j1 = Nt - l;
        Pl(j0:j1, li) = sum(hh(j0:j1,:) .* conj(vv(j0+l:j1+l,:)), 2);
      else
        j0 = 1 - l; j1 = Nt;
        Pl(j0:j1, li) = sum(hh(j0:j1,:) .* conj(vv(j0+l:j1+l,:)), 2);
      end
    end
    cPl = [complex(zeros(1, nlag)); cumsum(double(Pl), 1)];
    cEh = [0; cumsum(double(Ehh))];
    cEv = [0; cumsum(double(Evv))];
    clear Pl hh vv;

    % interpft resamples the nlag-point correlation over the SAME period, so
    % the fine grid starts at -maxlag with spacing 1/PAD and has nlag*PAD
    % points. Treating the correlation as periodic is safe here because it
    % has decayed to noise well inside +-8 m of lag - the column accumulates
    % ~3.8 m - so the wrap contributes nothing near the peak.
    nfine = nlag * PAD;
    lag_fine = -maxlag + (0:nfine-1).' / PAD;

    nkeep = nkeep + 1;
    L(nkeep).tag = tag;
    L(nkeep).lat = lat_n;
    L(nkeep).lon = lon_n;
    L(nkeep).nnear = numel(near);
    L(nkeep).t_surf = t_surf;
    L(nkeep).seg_list = SEG_LIST;

    for vi = 1:numel(SEG_LIST)
      SEG_M = SEG_LIST(vi);
      STEP_M = SEG_M / 4;      % their 12 m / 3 m ratio, held across variants
      zseg = (Z_TOP:STEP_M:Z_BOT-SEG_M).';
      nseg = numel(zseg);
      zc = zseg + SEG_M/2;
      dt_hat = nan(nseg,1); coh = nan(nseg,1);
      prev = NaN; z_prev = NaN;   % their surface-downward continuity track
      for si = 1:nseg
        i0 = find(z_raw >= zseg(si), 1);
        i1 = find(z_raw < zseg(si)+SEG_M, 1, 'last');
        if isempty(i0) || isempty(i1) || (i1-i0+1) < 16, continue; end
        if i0-maxlag < 1 || i1+maxlag > Nt, continue; end
        num = (cPl(i1+1,:) - cPl(i0,:)).';             % lag scan, complex
        eh = cEh(i1+1) - cEh(i0);
        ev = cEv(i1+1+lags) - cEv(i0+lags);            % VV window moves
        den = sqrt(eh * ev);
        c = num ./ max(den, realmin);

        % Interpolate the correlation to the x PAD lag grid - the analogue
        % of their padding factor p, applied where it is cheap.
        cf = interpft(c, nfine);
        [cbest, k] = max(abs(cf));
        if ~isfinite(cbest), continue; end
        dtau = lag_fine(k)*dt_s + angle(cf(k))/(2*pi*FC);
        % Re-anchor to the cycle nearest the running estimate: the peak is
        % picked on |c|, which is cycle-blind, so a hop appears as a whole
        % 1/fc jump. This is the discrete form of their downward track.
        % Only across a SHORT gap: dt drifts ~0.032 ns/m at dlam 0.5, so
        % beyond GAP_M the running value is more than half a cycle stale and
        % re-anchoring to it would lock in the wrong cycle for good.
        if isfinite(prev) && (zc(si) - z_prev) <= GAP_M
          dtau = dtau - round((dtau - prev)*FC)/FC;
        end
        dt_hat(si) = dtau; coh(si) = cbest;
        if cbest >= COH_MIN, prev = dtau; z_prev = zc(si); end
      end

      % --- their quality gate, then the 100 m smoothing
      dt_hat(~isfinite(coh) | coh < COH_MIN) = NaN;
      nsm = max(3, round(SMOOTH_M / STEP_M));
      dt_sm = movmean(dt_hat, nsm, 'omitnan');
      dt_sm(~isfinite(dt_hat)) = NaN;

      % --- 200 m moving-window slope, then dlam
      nde = max(3, round(DERIV_M / STEP_M));
      half = floor(nde/2);
      grad = nan(nseg,1);
      for si = 1:nseg
        lo = max(1, si-half); hi = min(nseg, si+half);
        zz = zc(lo:hi); yy = dt_sm(lo:hi);
        m = isfinite(yy);
        if nnz(m) < 5, continue; end
        zz = zz(m); yy = yy(m);
        if (max(zz) - min(zz)) < 0.5*DERIV_M, continue; end
        pf = polyfit(zz, yy, 1);
        grad(si) = pf(1);
      end
      dlam = grad / GRAD_PER_DLAM;

      fprintf(['  seg %4.1f m: %d/%d pass coh %.2f; median coh %.2f; ' ...
               'dlam %.3f..%.3f\n'], SEG_M, nnz(isfinite(dt_hat)), nseg, ...
        COH_MIN, median(coh(isfinite(coh))), min(dlam), max(dlam));

      L(nkeep).z{vi} = zc;
      L(nkeep).dt{vi} = dt_hat;
      L(nkeep).dt_sm{vi} = dt_sm;
      L(nkeep).coh{vi} = coh;
      L(nkeep).dlam{vi} = dlam;
    end
    clear cPl cEh cEv;
  catch ME
    % Warn and continue: one unusable frame must not cost the other eight,
    % and the azimuthal solve downstream needs only three azimuths.
    warning('egrip_zeising:frameFailed', '%s failed (%s): %s', ...
      tag, ME.identifier, ME.message);
    clear H V hh vv Pl cPl cEh cEv;
    continue;
  end
end

if nkeep < 3
  error('egrip_zeising:tooFewLines', ...
    ['only %d of %d lines produced a profile; the two-parameter azimuthal ' ...
     'fit needs at least 3'], nkeep, size(frames,1));
end

params = struct('SEG_LIST', SEG_LIST, 'COH_MIN', COH_MIN, ...
  'SMOOTH_M', SMOOTH_M, 'DERIV_M', DERIV_M, 'PAD', PAD, 'DEPS', DEPS, ...
  'FC', FC, 'RADIUS_KM', RADIUS_KM, 'GRAD_PER_DLAM', GRAD_PER_DLAM, ...
  'MAXLAG_M', MAXLAG_M, 'GAP_M', GAP_M);
out_fn = fullfile(out_dir, 'egrip_zeising.mat');
save(out_fn, '-v7.3', 'L', 'params');
fprintf('\nwrote %s (%d lines)\n', out_fn, nkeep);
