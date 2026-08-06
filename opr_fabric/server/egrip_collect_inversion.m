%EGRIP_COLLECT_INVERSION Gather our inversion's dlam for the EastGRIP lines.
%
% run_negis_fabric.m, driven with the nine borehole-proximal frames, leaves
% one CSARP_fabric_deltak_negis product per line. Each is small already
% (10 intervals x ~3 blocks), so this just concatenates the fields the
% method comparison needs into a single file to mirror back, rather than
% copying nine products.
%
% This is the path that carries AMPLITUDE: fabric_task was run with
% inversion = 'joint', i.e. ptt.invertHorizontalFabricJoint, which fits the
% returned power alongside dtau. That is the distinction being tested in
% scripts/figures/egrip_method_compare.py against the Zeising et al. (2023)
% phase-only estimator, so the two must come from the same lines, the same
% cull and the same surface pick - which they do, because both start from
% the frame list in egrip_zeising.m.
scratch = '/kucresis/scratch/hoffmana_sta/fabric';
season = '2024_Greenland_Ground2';
in_root = fullfile(scratch, season, 'CSARP_fabric_deltak_negis');
out_fn = fullfile(scratch, 'stages', 'egrip_inversion.mat');

frames = { '20240628_01', 1; '20240626_03', 5; '20240626_01', 1; ...
           '20240619_01', 1; '20240621_01', 1; '20240620_01', 1; ...
           '20240621_01', 10; '20240619_01', 5; '20240620_02', 3 };

R = struct([]);
n = 0;
for fi = 1:size(frames,1)
  day_seg = frames{fi,1};
  frm = frames{fi,2};
  tag = sprintf('%s_%03d', day_seg, frm);
  fn = fullfile(in_root, day_seg, sprintf('Data_%s_%03d.mat', day_seg, frm));
  if exist(fn, 'file') ~= 2
    warning('egrip_collect:missing', 'no inversion output for %s', tag);
    continue;
  end
  d = load(fn, 'dlam', 'dlam_top_depth', 'dlam_bot_depth', 'dlam_quality', ...
    'dlam_clipped', 'dlam_interpolated', 'Latitude', 'Longitude', ...
    'dtau_rms', 'phase_sign');
  n = n + 1;
  R(n).tag = tag;
  R(n).dlam = single(d.dlam);
  R(n).top = single(d.dlam_top_depth);
  R(n).bot = single(d.dlam_bot_depth);
  R(n).quality = single(d.dlam_quality);
  R(n).clipped = single(d.dlam_clipped);
  R(n).interpolated = single(d.dlam_interpolated);
  R(n).lat = d.Latitude;
  R(n).lon = d.Longitude;
  R(n).rms = d.dtau_rms;
  R(n).phase_sign = d.phase_sign;
  fprintf('%-18s %d intervals x %d blocks, %.0f-%.0f m, rms %.3f ns\n', ...
    tag, size(d.dlam,1), size(d.dlam,2), min(d.dlam_top_depth(:)), ...
    max(d.dlam_bot_depth(:)), median(d.dtau_rms));
end

if n < 3
  error('egrip_collect:tooFew', ...
    'only %d of %d lines had inversion output', n, size(frames,1));
end
save(out_fn, '-v7.3', 'R');
fprintf('\nwrote %s (%d lines)\n', out_fn, n);
