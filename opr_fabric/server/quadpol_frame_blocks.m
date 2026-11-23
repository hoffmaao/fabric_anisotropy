function [Dc, meta, z] = quadpol_frame_blocks(day_seg, frm, cache_dir, site_root, blk_m, z_max, nsub)
%QUADPOL_FRAME_BLOCKS One frame's coregistered channels, cut into along-track blocks.
%
% [Dc, meta, z] = quadpol_frame_blocks(day_seg, frm, cache_dir, site_root, blk_m, z_max, nsub)
%
% Reads the frame's coregistered channels (creg_<tag>.mat in cache_dir, as
% run_quadpol_pipeline caches them) down to z_max, rebuilds the per-trace
% heading the way the pipeline does (great-circle azimuth over SMH traces,
% smoothed on the doubled-angle phasor) from the season's polarimetric
% product positions, and cuts the frame into blocks of ~blk_m metres. Dc
% is a cell of channel structs (hh, vv, hv, vh, heading in rad), the form
% run_node_frame hands to nodeFabric; blocks of fewer than 2 nsub traces
% are left out. meta has per block x (along-track centre, m), lat, lon,
% az (mean heading, deg, mod 180), n (traces) and used (the block's index
% in the frame's full cut).
CH = {'hh', 'vv', 'hv', 'vh'};
day_seg = char(day_seg);
tag = sprintf('%s_%03d', day_seg, frm);
c = load(fullfile(cache_dir, ['creg_' tag '.mat']), 'hh', 'vv', 'hv', 'vh', 'z');
z = double(c.z(:));
keep = z <= z_max;
T = struct();
for k = 1:4, T.(CH{k}) = c.(CH{k})(keep, :); end
z = z(keep);
clear c
Nx = size(T.hh, 2);

name = sprintf('Data_%s.mat', tag);
pol_fn = fullfile(site_root, 'CSARP_polarimetric', day_seg, name);
if exist(pol_fn, 'file') ~= 2
  pol_fn = fullfile(site_root, 'CSARP_polarimetric_unwrap', day_seg, name);
end
Pq = load(pol_fn, 'Latitude', 'Longitude');
la = double(Pq.Latitude(:)); lo = double(Pq.Longitude(:));
if numel(la) ~= Nx
  error('quadpol_frame_blocks:positions', '%s: %d positions for %d traces', tag, numel(la), Nx);
end
R_E = 6371000;
dph = deg2rad(diff(la)); dlo = deg2rad(diff(lo));
aa = sin(dph/2).^2 + cos(deg2rad(la(1:end-1))) .* cos(deg2rad(la(2:end))) .* sin(dlo/2).^2;
x_along = [0; cumsum(2 * R_E * asin(min(1, sqrt(aa))))];
SMH = min(51, max(3, floor(Nx / 10)));
ih0 = 1:(Nx-SMH); ih1 = (1+SMH):Nx;
p0 = deg2rad(la(ih0)); p1 = deg2rad(la(ih1)); dl = deg2rad(lo(ih1) - lo(ih0));
azr = mod(rad2deg(atan2(sin(dl).*cos(p1), cos(p0).*sin(p1) - sin(p0).*cos(p1).*cos(dl))), 180);
xi = (1:numel(azr)).' + SMH/2;
ph = interp1(xi, exp(2i*deg2rad(azr)), min(max((1:Nx).', xi(1)), xi(end)), 'linear');
az_tr = mod(rad2deg(angle(ph))/2, 180);

nb = max(1, round(x_along(end) / blk_m));
edges = linspace(0, x_along(end), nb + 1);
meta = struct('x', [], 'lat', [], 'lon', [], 'az', [], 'n', [], 'used', []);
Dc = {};
for b = 1:nb
  js = find(x_along >= edges(b) & (x_along < edges(b+1) | b == nb));
  if numel(js) < 2 * nsub, continue; end
  meta.x(end+1) = mean(x_along(js));
  meta.lat(end+1) = mean(la(js)); meta.lon(end+1) = mean(lo(js));
  meta.az(end+1) = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(js))))))/2, 180);
  meta.n(end+1) = numel(js);
  meta.used(end+1) = b;
  S = struct();
  for k = 1:4, S.(CH{k}) = double(T.(CH{k})(:, js)); end
  S.heading = deg2rad(az_tr(js)).';
  Dc{end+1} = S; %#ok<AGROW>
end
end
