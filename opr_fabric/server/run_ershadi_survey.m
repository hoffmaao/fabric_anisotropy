%RUN_ERSHADI_SURVEY Ershadi et al. (2022) fabric inversion, per heading block.
%
% Applies ptt.ershadiFabric to every frame of a survey, in BLOCKS of
% near-constant heading rather than a frame at a time.
%
% WHY BLOCKS. Ershadi et al. process a stationary quad-pol sounding: one
% place, one antenna orientation, many chirps. A towed survey is not that.
% Averaging a whole frame in the antenna frame silently violates the
% assumption, and does so in a direction that matters - a signal fixed in
% the ICE rotates through the antenna frame as the vehicle wanders and
% averages down, while a signal fixed in the ANTENNAS adds coherently. On
% Ridge A 45 of 47 frames wander by more than 5 deg, with 5-95 spreads
% reaching 29, so the frame-at-a-time average was biased toward an
% antenna-fixed answer by construction. A 200-trace block is ~200 m of
% track, over which the measured wander is under a degree, which is the
% stationary case the published method assumes.
%
% THE HEADING TEST, at block level. Each block yields (heading, theta) and
% the blocks span the survey's full range of headings, so the question
% "does the recovered orientation follow the ice or the antennas?" is asked
% with hundreds of samples instead of tens, and asked WITHIN frames as well
% as between them - the same ice, the same calibration, the same
% processing, only the heading differing.
%
%   matlab -batch "site='ridge_a'; run_ershadi_survey"
if ~exist('site', 'var') || isempty(site), site = 'ridge_a'; end
switch site
  case 'ridge_a'
    site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
  case 'taylor_dome'
    site_root = '/cresis/dataproducts/opr_data/accum/2025_Antarctica_Ground2';
  case 'thwaites'
    site_root = '/cresis/dataproducts/opr_data/accum/2023_Antarctica_Ground';
  case 'eastwind'
    site_root = '/cresis/dataproducts/opr_data/accum/2022_Antarctica_Ground';
  otherwise
    error('run_ershadi_survey:site', 'unknown site %s', site);
end
% Repo root and <work> are derived from this script's own location - see
% fabric_paths.
[code, scratch] = fabric_paths();
out_fn = fullfile(scratch, 'stages', sprintf('ershadi_%s.mat', site));
addpath(code);

CHAN = {'hh','vv','hv','vh'};
NBLK = 200;             % traces per heading block
FC = 750e6;
C0 = 299792458;
EPS_ICE = 3.171;
C_ICE = C0 / sqrt(EPS_ICE);
Z_MAX = 1500;
Z_BAND = [200 1200];    % depth band the per-block summary is taken over
SMH = 51;               % baseline for the per-trace heading

hh_dir = fullfile(site_root, 'CSARP_standardphase_HH');
d = dir(fullfile(hh_dir, '*', 'Data_*.mat'));
d = d(cellfun(@isempty, regexp({d.name}, '^Data_img_', 'once')));
if isempty(d)
  error('run_ershadi_survey:noFrames', 'no frames under %s', hh_dir);
end
fprintf('site %s: %d frames\n', site, numel(d));

B = struct([]);
nb = 0;
for fi = 1:numel(d)
  tag = regexprep(d(fi).name, '^Data_|\.mat$', '');
  [~, ds] = fileparts(d(fi).folder);
  fprintf('\n=== %s (%d of %d) ===\n', tag, fi, numel(d));
  try
    S = struct(); Tv = []; ref = struct();
    for k = 1:4
      fn = fullfile(site_root, ['CSARP_standardphase_' upper(CHAN{k})], ...
        ds, d(fi).name);
      if exist(fn, 'file') ~= 2
        error('run_ershadi_survey:missing', 'missing %s channel', CHAN{k});
      end
      q = load(fn, 'Data', 'Time', 'Latitude', 'Longitude', 'Surface');
      if isreal(q.Data)
        error('run_ershadi_survey:real', '%s is real-valued', CHAN{k});
      end
      if k == 1
        Tv = q.Time(:); ref = q;
      elseif ~isequal(size(q.Data), size(ref.Data)) ...
          || max(abs(q.Time(:) - Tv)) > 1e-12
        error('run_ershadi_survey:axes', '%s axes differ from HH', CHAN{k});
      end
      S.(CHAN{k}) = q.Data;
    end
    [~, Nx] = size(S.hh);

    surf_t = median(ref.Surface(:), 'omitnan');
    if ~isfinite(surf_t)
      pw = abs(S.hh).^2;
      thr = 0.05 * max(pw, [], 1);
      st = nan(1, Nx);
      for j = 1:Nx
        i0 = find(pw(:,j) > thr(j), 1);
        if ~isempty(i0), st(j) = Tv(i0); end
      end
      surf_t = median(st, 'omitnan');
      clear pw thr;
    end
    z = (Tv - surf_t) * C_ICE / 2;
    keep = z >= 0 & z <= Z_MAX;
    for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(keep, :); end
    z = z(keep);
    zb = z > Z_BAND(1) & z < Z_BAND(2);

    la = ref.Latitude(:); lo = ref.Longitude(:);
    ih0 = 1:(Nx-SMH); ih1 = (1+SMH):Nx;
    ph0 = deg2rad(la(ih0)); ph1 = deg2rad(la(ih1));
    dlh = deg2rad(lo(ih1) - lo(ih0));
    az_tr = mod(rad2deg(atan2(sin(dlh).*cos(ph1), ...
      cos(ph0).*sin(ph1) - sin(ph0).*cos(ph1).*cos(dlh))), 180);
    az_tr = interp1((1:numel(az_tr)).' + SMH/2, az_tr, (1:Nx).', ...
      'linear', 'extrap');

    nblocks = max(1, floor(Nx / NBLK));
    kept = 0;
    for b = 1:nblocks
      j0 = (b-1)*NBLK + 1;
      j1 = min(b*NBLK, Nx);
      if j1 - j0 < 32, continue; end
      Sb = struct();
      for k = 1:4, Sb.(CHAN{k}) = S.(CHAN{k})(:, j0:j1); end
      % heading of this block, and how constant it is across it
      ab = az_tr(j0:j1);
      az_blk = mod(rad2deg(angle(mean(exp(2i*deg2rad(ab)))))/2, 180);
      blk_spread = diff(prctile(mod(ab - az_blk + 90, 180) - 90, [5 95]));

      o = ptt.ershadiFabric(Sb, z, struct('fc', FC, 'psi_step_deg', 1, ...
        'win_m', 30, 'grad_win_m', 25, 'coh_min', 0.4, 'deramped', true));

      th = rad2deg(o.theta(zb));
      good = isfinite(th);
      if nnz(good) < 20, clear Sb; continue; end
      th_ant = mod(rad2deg(angle(mean(exp(2i*deg2rad(th(good))))))/2, 180);
      pv = mean(abs(Sb.hv(:)).^2) / max(mean(abs(Sb.hh(:)).^2), realmin);
      rc = abs(mean(Sb.hv(:) .* conj(Sb.vh(:)))) / ...
        sqrt(max(mean(abs(Sb.hv(:)).^2)*mean(abs(Sb.vh(:)).^2), realmin));

      nb = nb + 1; kept = kept + 1;
      B(nb).tag = tag;
      B(nb).blk = b;
      B(nb).az = az_blk;
      B(nb).az_spread = blk_spread;
      B(nb).theta_ant = th_ant;
      B(nb).theta_geo = mod(th_ant + az_blk, 180);
      B(nb).dlam = median(o.dlam(zb), 'omitnan');
      B(nb).cmag = median(o.Cmag(zb, 1), 'omitnan');
      B(nb).qual = median(o.theta_quality(zb), 'omitnan');
      B(nb).xr = 10*log10(max(pv, realmin));
      B(nb).recip = rc;
      B(nb).lat = median(la(j0:j1)); B(nb).lon = median(lo(j0:j1));
      clear Sb o;
    end
    fprintf('  %d blocks kept (heading %.0f deg)\n', kept, ...
      mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr)))))/2, 180));
    save(out_fn, '-v7.3', 'B');
    clear S;
  catch ME
    warning('run_ershadi_survey:frameFailed', '%s failed (%s): %s', ...
      tag, ME.identifier, ME.message);
    clear S Sb;
    continue;
  end
end

if nb < 1
  error('run_ershadi_survey:tooFew', 'no blocks produced');
end
fprintf('\nwrote %s (%d blocks from %d frames)\n', out_fn, nb, numel(d));
