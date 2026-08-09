%RUN_QUADPOL_FRAME Quad-pol fabric inversion on one real frame.
%
% Loads all four channels of a CSARP_standardphase_{HH,VV,HV,VH} frame,
% forms the scattering-matrix moments, synthesizes the azimuth sweep and
% inverts for orientation and contrast (ptt.quadpolMoments ->
% quadpolAzimuth -> quadpolFabric).
%
% CALIBRATION IS CHECKED BEFORE THE SCIENCE, and is reported whether or not
% it passes. The whole method rests on the cross-polarized nulls being
% real; finite antenna isolation leaks co-pol into the HV/VH channels and
% fills those nulls in, which would bias the recovered orientation toward
% whatever the leakage pattern is and flatten the modulation depth that
% the quality gate keys on. Two diagnostics:
%   1. RECIPROCITY. A reciprocal medium has S_hv == S_vh, so any measured
%      difference is system, not ice. Reported as the coherence between
%      the two channels and their power ratio.
%   2. CROSS-POL RATIO. |S_hv|^2 / |S_hh|^2. If this sits at a constant
%      floor with depth rather than oscillating, it is leakage rather than
%      birefringence, because real cross-pol must go through nulls where
%      delta passes multiples of 2*pi.
%
% Set `day_seg`, `frm` and `site_root` before running, or take the
% defaults (Ridge A leg 20250108_02_009, the frame run_sections.m sections
% and the SCAR figure draws, so the quad-pol answer can be set against a
% result that already exists).
if ~exist('site_root', 'var') || isempty(site_root)
  site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
end
if ~exist('day_seg', 'var') || isempty(day_seg), day_seg = '20250108_02'; end
if ~exist('frm', 'var') || isempty(frm), frm = 9; end
if ~exist('out_fn', 'var') || isempty(out_fn)
  out_fn = sprintf('/kucresis/scratch/hoffmana_sta/fabric/stages/quadpol_%s_%03d.mat', ...
    day_seg, frm);
end
addpath('/kucresis/scratch/hoffmana_sta/fabric/code');

CHAN = {'hh','vv','hv','vh'};
NR = 9;                 % range looks for the moment matrix
PSI_STEP_DEG = 2;       % azimuth sweep spacing
FC = 750e6;
C0 = 299792458;
EPS_ICE = 3.171;        % ptt.constants eps_bar
C_ICE = C0 / sqrt(EPS_ICE);
Z_MAX = 1500;

name = sprintf('Data_%s_%03d.mat', day_seg, frm);
S = struct(); Tv = []; ref = struct();
for k = 1:4
  fn = fullfile(site_root, ['CSARP_standardphase_' upper(CHAN{k})], ...
    day_seg, name);
  if exist(fn, 'file') ~= 2
    error('run_quadpol_frame:missing', 'missing %s', fn);
  end
  q = load(fn, 'Data', 'Time', 'GPS_time', 'Latitude', 'Longitude', ...
    'Surface', 'Elevation');
  if isreal(q.Data)
    error('run_quadpol_frame:real', ...
      '%s is real-valued; the scattering matrix needs phase', fn);
  end
  if k == 1
    Tv = q.Time(:);
    ref = q;
  else
    if ~isequal(size(q.Data), size(ref.Data))
      error('run_quadpol_frame:size', '%s size differs from HH', CHAN{k});
    end
    if max(abs(q.Time(:) - Tv)) > 1e-12
      error('run_quadpol_frame:time', '%s range axis differs from HH', CHAN{k});
    end
  end
  S.(CHAN{k}) = q.Data;
end
[Nt, Nx] = size(S.hh);
fprintf('%s_%03d: %d samples x %d traces\n', day_seg, frm, Nt, Nx);

% --- depth axis from the product's own surface pick
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
  fprintf('  Surface was NaN; picked leading edge at %.3f us\n', surf_t*1e6);
  clear pw thr;
end
z = (Tv - surf_t) * C_ICE / 2;
keep = z >= 0 & z <= Z_MAX;
for k = 1:4, S.(CHAN{k}) = S.(CHAN{k})(keep, :); end
z = z(keep);
fprintf('  depth window %.0f..%.0f m (%d samples, dz %.2f m)\n', ...
  z(1), z(end), numel(z), median(diff(z)));

% --- calibration diagnostics, before anything is inverted
p_hh = mean(abs(S.hh).^2, 2);
p_hv = mean(abs(S.hv).^2, 2);
p_vh = mean(abs(S.vh).^2, 2);
c_x = mean(S.hv .* conj(S.vh), 2) ./ ...
  sqrt(max(p_hv .* p_vh, realmin));
mid = z > 200 & z < 1200;
fprintf('\n  --- calibration ---\n');
fprintf('  HV/VH coherence      median %.3f  (1 = reciprocal)\n', ...
  median(abs(c_x(mid)), 'omitnan'));
fprintf('  HV/VH power ratio    median %.2f dB (0 = balanced)\n', ...
  10*log10(median(p_hv(mid)./max(p_vh(mid), realmin), 'omitnan')));
xr = 10*log10(0.5*(p_hv + p_vh) ./ max(p_hh, realmin));
fprintf('  cross/co power       median %.1f dB, 10-90%% spread %.1f dB\n', ...
  median(xr(mid), 'omitnan'), diff(prctile(xr(mid), [10 90])));
fprintf('  (a flat cross/co with depth is leakage; real cross-pol must\n');
fprintf('   null wherever delta passes a multiple of 2*pi)\n\n');

% --- moments, sweep, inversion
M = ptt.quadpolMoments(S, [NR Nx]);
psi = (0:PSI_STEP_DEG:180-PSI_STEP_DEG) * pi/180;
A = ptt.quadpolAzimuth(M, psi);
out = ptt.quadpolFabric(A, z, struct('fc', FC, 'win_m', 50, ...
  'grad_win_m', 200));

% --- antenna frame -> geographic. H is along-track on this system, so the
% recovered angle is relative to the heading and the track azimuth has to
% be added back before it can be compared with anything geographic.
la = ref.Latitude(:); lo = ref.Longitude(:);
p0 = deg2rad(la(1)); p1 = deg2rad(la(end));
dl = deg2rad(lo(end) - lo(1));
track_az = mod(rad2deg(atan2(sin(dl)*cos(p1), ...
  cos(p0)*sin(p1) - sin(p0)*cos(p1)*cos(dl))), 180);
theta_geo = mod(rad2deg(out.theta) + track_az, 180);

fprintf('  track azimuth        %.1f deg E of N\n', track_az);
fprintf('  aniso (modulation)   median %.3f\n', median(out.aniso(mid), 'omitnan'));
fprintf('  theta (antenna)      median %.1f deg\n', ...
  mod(rad2deg(median(out.theta(mid), 'omitnan')), 90));
fprintf('  theta (geographic)   median %.1f deg E of N\n', ...
  mod(median(theta_geo(mid), 'omitnan'), 180));
fprintf('  dlam  (phase grad)   median %.3f\n', median(out.dlam(mid), 'omitnan'));
fprintf('  dlam  (cross nodes)  median %.3f\n', ...
  median(out.dlam_node(mid), 'omitnan'));
fprintf('  branch_flipped       %d\n', out.branch_flipped);

save(out_fn, '-v7.3', 'out', 'z', 'psi', 'A', 'theta_geo', 'track_az', ...
  'day_seg', 'frm', 'c_x', 'xr', 'p_hh', 'p_hv', 'p_vh');
fprintf('\nwrote %s\n', out_fn);
