%NEGIS_FULL_INVENTORY Every frame of the season, and whether it is usable.
%
% "Usable" means the qlook product is COMPLEX (written with inc_dec = 0, so
% the HH/VV phase difference survives) and its day GPS file covers the
% frame. Both are per-frame facts and neither is recorded anywhere, so the
% full-season run needs them measured before it is queued rather than
% discovered 11 hours in.
season = '/cresis/dataproducts/opr_data/accum/2024_Greenland_Ground2';
gps_dir = '/cresis/dataproducts/opr_data/opr_support/gps/2024_Greenland_Ground2';
hh_root = fullfile(season, 'CSARP_qlook_HH');
vv_root = fullfile(season, 'CSARP_qlook_VV');

segs = dir(hh_root);
segs = segs([segs.isdir] & ~startsWith({segs.name}, '.'));
fprintf('%d segments\n', numel(segs));

tot = 0; ok_n = 0;
lines = {};
for si = 1:numel(segs)
  ds = segs(si).name;
  d = dir(fullfile(hh_root, ds, 'Data_*.mat'));
  tot = tot + numel(d);
  % inc_dec is a per-segment processing choice, so the first frame settles
  % whether the whole segment carries phase - no need to open all of them.
  if isempty(d), continue; end
  try
    q = load(fullfile(d(1).folder, d(1).name), 'Data');
    cpx = ~isreal(q.Data);
  catch ME
    fprintf('  %-14s UNREADABLE (%s)\n', ds, ME.message);
    continue;
  end
  clear q;
  gfn = fullfile(gps_dir, sprintf('gps_%s.mat', ds(1:8)));
  has_gps = exist(gfn, 'file') == 2;
  nvv = numel(dir(fullfile(vv_root, ds, 'Data_*.mat')));
  fprintf('  %-14s %3d HH / %3d VV frames  complex=%d  gps=%d\n', ds, ...
    numel(d), nvv, cpx, has_gps);
  if ~cpx || ~has_gps, continue; end
  for k = 1:numel(d)
    t = regexp(d(k).name, 'Data_(\d{8}_\d{2})_(\d{3})\.mat', 'tokens', 'once');
    if isempty(t), continue; end
    if exist(fullfile(vv_root, ds, d(k).name), 'file') ~= 2, continue; end
    ok_n = ok_n + 1;
    lines{end+1} = sprintf('''%s'',%d', t{1}, str2double(t{2})); %#ok<SAGROW>
  end
end
fprintf('\n%d frames total, %d complex with VV and GPS\n', tot, ok_n);
fprintf('\ntargets = { ...\n');
for k = 1:8:numel(lines)
  fprintf(' %s; ...\n', strjoin(lines(k:min(k+7,numel(lines))), ';'));
end
fprintf('};\n');
