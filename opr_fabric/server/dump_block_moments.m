%DUMP_BLOCK_MOMENTS Along-track block moment matrices of a frame, from its coreg cache.
%
%   matlab -batch "day_seg='20240120_06'; frm=6; blk_m=1000; dump_block_moments"
%
% Writes stages/quadpol/moments/blkmom_<tag>_B<blk_m>.mat holding, per
% along-track block of ~blk_m metres, the 4x4 polarimetric moment matrix
% of the coregistered channels [Nt x 4 x 4 x nb] (single), in BOTH the
% antenna frame (as recorded) and the geographic frame (heading-rotated
% per 200-trace sub-block, as ptt.quadpolFrameTheta pools them), with the
% block centre position, heading and trace count.
%
% WHY. Every estimator in +ptt that works on moments - the coherence LS
% (ptt.quadpolFabricLS), the power-extinction fit
% (ptt.quadpolFabricPower), and whatever comes next - can then run on a
% frame in seconds on a laptop instead of a 30-70 min pipeline pass on
% the shared node, and on exactly the same blocks, which is what a
% comparison between them needs. A block matrix is ~1.3 MB, so a frame is
% tens of MB against the 0.7 GB cache.
%
% POSITIONS come from the frame's shipped polarimetric product where the
% season has one, and otherwise (EastGRIP qlook) from the existing
% quadpol_section product's per-block sec_lat/sec_lon/sec_az, so the dump
% never repeats the pipeline's trajectory rebuild and never disagrees with
% it. Along-track distance is haversine-cumulative, as in the pipeline.
%
% Requires: the coreg cache creg_<tag>.mat (run the pipeline first).
[fabric_code, fabric_work] = fabric_paths();
addpath(fabric_code);
if ~exist('day_seg', 'var'), error('set day_seg'); end
if ~exist('frm', 'var'), error('set frm'); end
day_seg = char(day_seg);
if ~exist('blk_m', 'var') || isempty(blk_m), blk_m = 1000; end
if ~exist('nblk_rot', 'var') || isempty(nblk_rot), nblk_rot = 200; end
tag = sprintf('%s_%03d', day_seg, frm);
stage = fullfile(fabric_work, 'stages', 'quadpol');
cache_fn = fullfile(stage, 'coreg_cache', sprintf('creg_%s.mat', tag));
prod_fn = fullfile(stage, sprintf('quadpol_section_%s.mat', tag));
out_dir = fullfile(stage, 'moments');
if exist(out_dir, 'dir') ~= 7, mkdir(out_dir); end
out_fn = fullfile(out_dir, sprintf('blkmom_%s_B%d.mat', tag, round(blk_m)));
CHAN = {'hh', 'vv', 'hv', 'vh'};

c = load(cache_fn);
T = struct('hh', c.hh, 'vv', c.vv, 'hv', c.hv, 'vh', c.vh);
z = double(c.z(:));
Nx = size(T.hh, 2);
clear c;

% --- per-trace positions and heading
la = []; lo = [];
if exist('site_root', 'var') && ~isempty(site_root)
  site_root = char(site_root);
  name = sprintf('Data_%s_%03d.mat', day_seg, frm);
  pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
  if exist(pol_fn, 'file') ~= 2
    pol_fn = fullfile(site_root, 'CSARP_polarimetric_unwrap', day_seg, name);
  end
  if exist(pol_fn, 'file') == 2
    P = load(pol_fn, 'Latitude', 'Longitude');
    la = double(P.Latitude(:)); lo = double(P.Longitude(:));
  end
end
pos_source = 'polarimetric product';
if numel(la) ~= Nx
  % fall back to the section product's block positions, interpolated to
  % traces: the qlook seasons have no polarimetric product, and a culled
  % trace axis would not match it anyway
  if exist(prod_fn, 'file') ~= 2
    error('dump_block_moments:positions', ...
      'no per-trace positions for %s: no polarimetric product and no section product %s', tag, prod_fn);
  end
  r = load(prod_fn, 'res'); r = r.res;
  nbt = double(r.nblk_tr);
  xb = ((1:numel(r.sec_lat)) - 0.5) * nbt;     % block centres in traces
  la = interp1(xb, double(r.sec_lat(:)), (1:Nx).', 'linear', 'extrap');
  lo = interp1(xb, double(r.sec_lon(:)), (1:Nx).', 'linear', 'extrap');
  pos_source = sprintf('section product blocks (%d traces)', nbt);
end
R_E = 6371000;
dph = deg2rad(diff(la)); dlo = deg2rad(diff(lo));
aa = sin(dph/2).^2 + cos(deg2rad(la(1:end-1))) .* cos(deg2rad(la(2:end))) .* sin(dlo/2).^2;
x_along = [0; cumsum(2 * R_E * asin(min(1, sqrt(aa))))];
% per-trace heading axis (deg mod 180), smoothed on the doubled-angle
% phasor as the pipeline does
SMH = min(51, max(3, floor(Nx / 10)));
ih0 = 1:(Nx-SMH); ih1 = (1+SMH):Nx;
p0 = deg2rad(la(ih0)); p1 = deg2rad(la(ih1)); dl = deg2rad(lo(ih1) - lo(ih0));
azr = mod(rad2deg(atan2(sin(dl).*cos(p1), cos(p0).*sin(p1) - sin(p0).*cos(p1).*cos(dl))), 180);
xi = (1:numel(azr)).' + SMH/2;
ph = interp1(xi, exp(2i*deg2rad(azr)), min(max((1:Nx).', xi(1)), xi(end)), 'linear');
az_tr = mod(rad2deg(angle(ph))/2, 180);

% --- blocks
nb = max(1, round(x_along(end) / blk_m));
edges = linspace(0, x_along(end), nb + 1);
Nt = numel(z);
M_ant = zeros(Nt, 4, 4, nb, 'single');
M_geo = zeros(Nt, 4, 4, nb, 'single');
blk = struct('x', nan(1, nb), 'lat', nan(1, nb), 'lon', nan(1, nb), ...
  'az', nan(1, nb), 'n', zeros(1, nb));
for b = 1:nb
  js = find(x_along >= edges(b) & (x_along < edges(b+1) | b == nb));
  blk.n(b) = numel(js);
  if numel(js) < 32, continue; end
  blk.x(b) = mean(x_along(js)); blk.lat(b) = mean(la(js)); blk.lon(b) = mean(lo(js));
  blk.az(b) = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(js))))))/2, 180);
  Sb = struct();
  for k = 1:4, Sb.(CHAN{k}) = T.(CHAN{k})(:, js); end
  M_ant(:, :, :, b) = single(ptt.quadpolMoments(Sb, [1 numel(js)]));
  Mg = zeros(Nt, 4, 4); wsum = 0;
  for jr = 1:nblk_rot:numel(js)
    jr1 = min(jr + nblk_rot - 1, numel(js));
    if jr1 - jr < 16, continue; end
    jj = js(jr:jr1);
    Sr = struct();
    for k = 1:4, Sr.(CHAN{k}) = T.(CHAN{k})(:, jj); end
    hb = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(jj))))))/2, 180);
    Mg = Mg + ptt.rotateMoments(ptt.quadpolMoments(Sr, [1 numel(jj)]), -deg2rad(hb)) * numel(jj);
    wsum = wsum + numel(jj);
  end
  if wsum > 0, M_geo(:, :, :, b) = single(Mg / wsum); end
  fprintf('block %2d/%d: %4d traces, x %6.0f m, heading %5.1f\n', b, nb, numel(js), blk.x(b), blk.az(b));
end
save(out_fn, '-v7.3', 'tag', 'z', 'M_ant', 'M_geo', 'blk', 'blk_m', 'nblk_rot', ...
  'pos_source', 'x_along', 'az_tr');
fprintf('wrote %s (%s)\n', out_fn, pos_source);
