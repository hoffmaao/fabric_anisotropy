%RUN_QUADPOL_SURVEY Quad-pol fabric inversion over every frame of a survey.
%
% Runs the ptt.quadpolMoments -> quadpolAzimuth -> quadpolFabric chain on
% all four standardphase channels of every frame, and stores the per-frame
% profiles plus the calibration diagnostics.
%
% THE POINT IS THE HEADING TEST. The claim being made for quad-pol is that
% orientation is recoverable from a SINGLE line, with no azimuth diversity
% required. If that is true, theta expressed GEOGRAPHICALLY must come out
% the same on every frame regardless of which way the vehicle was driving;
% if instead the recovered angle tracks the track heading, the method is
% reporting the antenna frame back to itself and has measured nothing. A
% raster survey drives its legs at two headings roughly 90 deg apart, which
% is exactly the control needed, and it costs no extra acquisition.
%
% The calibration diagnostics are recorded per frame for the second open
% question: whether the cross-polarized pedestal that defeated the
% fringe-rate estimator on 20250108_02_009 is a fixed instrument property
% (antenna isolation - constant across frames, and removable) or a
% property of the ice (volume scattering - varying with site and depth).
% One frame cannot separate those; a survey can.
%
% Pick the survey with `site`, default ridge_a:
%   matlab -batch "site='ridge_a'; run_quadpol_survey"
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
    error('run_quadpol_survey:site', 'unknown site %s', site);
end
% Repo root and <work> are derived from this script's own location - see
% fabric_paths.
[code, scratch] = fabric_paths();
out_fn = fullfile(scratch, 'stages', sprintf('quadpol_%s.mat', site));
addpath(code);

CHAN = {'hh','vv','hv','vh'};
NR = 9;                 % range looks for the moment matrix
PSI_STEP_DEG = 2;
FC = 750e6;
C0 = 299792458;
EPS_ICE = 3.171;
C_ICE = C0 / sqrt(EPS_ICE);
Z_MAX = 1500;
ZSUB = 4;               % store every 4th depth sample; ~1.1 m is far finer
                        % than the 200 m gradient window resolves anyway
% Traces per heading block. The moments are averaged within a block in the
% ANTENNA frame and the block is then rotated to geographic before being
% accumulated, so the block only has to be short enough that the heading is
% constant across it. At ~1 m trace spacing 200 traces is ~200 m of track,
% over which the measured wander is under a degree.
NBLK = 200;

hh_dir = fullfile(site_root, 'CSARP_standardphase_HH');
d = dir(fullfile(hh_dir, '*', 'Data_*.mat'));
% Drop the img_NN waveform sub-images: they are the same geometry as their
% parent frame, so including them would triple the run and triple-count
% every frame in any statistic computed downstream.
d = d(cellfun(@isempty, regexp({d.name}, '^Data_img_', 'once')));
if isempty(d)
  error('run_quadpol_survey:noFrames', 'no frames under %s', hh_dir);
end
fprintf('site %s: %d frames under %s\n', site, numel(d), hh_dir);

Q = struct([]);
n = 0;
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
        error('run_quadpol_survey:missing', 'missing %s channel', CHAN{k});
      end
      q = load(fn, 'Data', 'Time', 'GPS_time', 'Latitude', 'Longitude', ...
        'Surface');
      if isreal(q.Data)
        error('run_quadpol_survey:real', '%s channel is real-valued', CHAN{k});
      end
      if k == 1
        Tv = q.Time(:); ref = q;
      elseif ~isequal(size(q.Data), size(ref.Data)) ...
          || max(abs(q.Time(:) - Tv)) > 1e-12
        error('run_quadpol_survey:axes', '%s axes differ from HH', CHAN{k});
      end
      S.(CHAN{k}) = q.Data;
    end
    [Nt, Nx] = size(S.hh);

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
    fprintf('  %d x %d, depth %.0f..%.0f m\n', Nt, Nx, z(1), z(end));

    % --- calibration, per frame
    p_hh = mean(abs(S.hh).^2, 2);
    p_hv = mean(abs(S.hv).^2, 2);
    p_vh = mean(abs(S.vh).^2, 2);
    recip = abs(mean(S.hv .* conj(S.vh), 2)) ./ sqrt(max(p_hv.*p_vh, realmin));
    xr = 10*log10(0.5*(p_hv + p_vh) ./ max(p_hh, realmin));
    mid = z > 200 & z < 1200;

    % --- per-trace heading, for the geographic-frame average
    la = ref.Latitude(:); lo = ref.Longitude(:);
    SMH = 51;
    ih0 = 1:(Nx-SMH); ih1 = (1+SMH):Nx;
    ph0 = deg2rad(la(ih0)); ph1 = deg2rad(la(ih1));
    dlh = deg2rad(lo(ih1) - lo(ih0));
    az_tr = mod(rad2deg(atan2(sin(dlh).*cos(ph1), ...
      cos(ph0).*sin(ph1) - sin(ph0).*cos(ph1).*cos(dlh))), 180);
    az_tr = interp1((1:numel(az_tr)).' + SMH/2, az_tr, (1:Nx).', ...
      'linear', 'extrap');

    % --- moments, both ways.
    % ANTENNA frame: average every trace as measured. This is what the
    % first run did and it prefers whatever is fixed to the instrument.
    % GEOGRAPHIC frame: average blocks of near-constant heading after
    % rotating each into a common north-referenced frame, which prefers
    % whatever is fixed in the ice. Comparing them is the only way to tell
    % an instrument leakage term from an artifact of the averaging.
    M = ptt.quadpolMoments(S, [NR Nx]);
    nb = max(1, floor(Nx / NBLK));
    Mg = complex(zeros(size(M)));
    wsum = 0;
    for b = 1:nb
      j0 = (b-1)*NBLK + 1;
      j1 = min(b*NBLK, Nx);
      if j1 - j0 < 8, continue; end
      Sb = struct();
      for k = 1:4, Sb.(CHAN{k}) = S.(CHAN{k})(:, j0:j1); end
      Mb = ptt.quadpolMoments(Sb, [NR (j1-j0+1)]);
      % rotate the block into geographic: antennas sit at az_blk, so the
      % rotation that carries them to north is -az_blk
      az_blk = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(j0:j1)))))) / 2, 180);
      Mg = Mg + ptt.rotateMoments(Mb, -deg2rad(az_blk)) * (j1-j0+1);
      wsum = wsum + (j1-j0+1);
      clear Sb Mb;
    end
    if wsum > 0, Mg = Mg / wsum; end
    hdg_spread = diff(prctile(mod(az_tr - median(az_tr) + 90, 180) - 90, ...
      [5 95]));
    clear S;
    psi = (0:PSI_STEP_DEG:180-PSI_STEP_DEG) * pi/180;
    A = ptt.quadpolAzimuth(M, psi);
    out = ptt.quadpolFabric(A, z, struct('fc', FC, 'win_m', 50, ...
      'grad_win_m', 200));
    Ag = ptt.quadpolAzimuth(Mg, psi);
    outg = ptt.quadpolFabric(Ag, z, struct('fc', FC, 'win_m', 50, ...
      'grad_win_m', 200));

    p0 = deg2rad(la(1)); p1 = deg2rad(la(end));
    dl = deg2rad(lo(end) - lo(1));
    track_az = mod(rad2deg(atan2(sin(dl)*cos(p1), ...
      cos(p0)*sin(p1) - sin(p0)*cos(p1)*cos(dl))), 180);
    theta_geo = mod(rad2deg(out.theta) + track_az, 180);
    % Already north-referenced: the block rotation removed the heading, so
    % adding track_az again would put it back.
    theta_geo_g = mod(rad2deg(outg.theta), 180);

    fprintf(['  track %5.1f (spread %4.1f) | theta_geo ant %5.1f  geo %5.1f' ...
             ' | dlam %.3f / %.3f | aniso %.2f / %.2f | x/co %.1f | ' ...
             'recip %.3f\n'], track_az, hdg_spread, ...
      mod(H_cmed(theta_geo(mid)), 180), mod(H_cmed(theta_geo_g(mid)), 180), ...
      median(out.dlam(mid), 'omitnan'), median(outg.dlam(mid), 'omitnan'), ...
      median(out.aniso(mid), 'omitnan'), median(outg.aniso(mid), 'omitnan'), ...
      median(xr(mid), 'omitnan'), median(recip(mid), 'omitnan'));

    s = 1:ZSUB:numel(z);
    n = n + 1;
    Q(n).tag = tag;
    Q(n).z = single(z(s));
    Q(n).theta_geo = single(theta_geo(s));
    Q(n).theta_ant = single(rad2deg(out.theta(s)));
    Q(n).dlam = single(out.dlam(s));
    Q(n).dlam_node = single(out.dlam_node(s));
    Q(n).aniso = single(out.aniso(s));
    Q(n).xr = single(xr(s));
    Q(n).recip = single(recip(s));
    Q(n).track_az = track_az;
    Q(n).lat = median(la, 'omitnan');
    Q(n).lon = median(lo, 'omitnan');
    Q(n).lat0 = la(1); Q(n).lon0 = lo(1);
    Q(n).lat1 = la(end); Q(n).lon1 = lo(end);
    Q(n).branch_flipped = out.branch_flipped;
    Q(n).hdg_spread = hdg_spread;
    Q(n).theta_geo_g = single(theta_geo_g(s));
    Q(n).dlam_g = single(outg.dlam(s));
    Q(n).aniso_g = single(outg.aniso(s));
    save(out_fn, '-v7.3', 'Q');   % incremental, as run_survey_fabric does
    clear M Mg A Ag out outg;
  catch ME
    warning('run_quadpol_survey:frameFailed', '%s failed (%s): %s', ...
      tag, ME.identifier, ME.message);
    clear S M Mg A Ag out outg;
    continue;
  end
end

if n < 1
  error('run_quadpol_survey:tooFew', 'no frames inverted');
end
fprintf('\nwrote %s (%d of %d frames)\n', out_fn, n, numel(d));

function m = H_cmed(a)
% Median of an axis-valued angle in degrees, modulo 180: take it on the
% doubled angle so 179 and 1 average to 0 rather than to 90.
a = a(isfinite(a));
if isempty(a), m = NaN; return; end
m = mod(0.5 * rad2deg(angle(mean(exp(1i*2*deg2rad(a))))), 180);
end
