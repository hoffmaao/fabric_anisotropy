%RUN_DELTAK_SCRATCH Delta-k (split-spectrum) fabric batch on the servers.
%   Runs fabric_task with dtau_source = 'deltak' (ptt.deltakTraveltime:
%   absolute dtau from the ref/sec SLC spectra, no phase unwrapping, no
%   fringe blending) over the products that carry the complex SLCs:
%     - 2024_Antarctica_Ground2 Paden CSARP_polarimetric (Ridge A grid):
%       the clean-data validation - the SNAPHU chain is healthy there
%       (dtau_rms ~0.2 ns, blend fringes 0..1, phase_sign -1), so
%       delta-k should reproduce fabric_joint_jp within noise while
%       confirming the sign and the ladder margins.
%     - 2023_Antarctica_Ground CSARP_polarimetric_unwrap (Thwaites line):
%       the margin frames where delta-k disagrees with the production
%       blend corrections and flips the sign.
%   Outputs land beside the existing fabric_joint/joint_jp results as
%   CSARP_fabric_deltak(_jp) for direct comparison (same block/interval
%   structure).
%
%   NOTE ptt.deltakTraveltime needs the UNREGISTERED sec: coregistration
%   shifts the envelope and removes the group delay. fabric_task loads
%   'ref' and 'sec' (never sec_reg) in deltak mode; products that saved
%   only sec_reg cannot be used.
%
%   Rerun-safe: frames whose output file already exists are skipped.
%
%   Launch on mem1 with:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('<code>/opr_fabric/server/run_deltak_scratch.m')"

% Repo root and <work> are derived from this script's own location - see
% fabric_paths.
[proj_root, scratch] = fabric_paths();
% {season, data_root, product name, out_path}
product_tbl = { ...
  '2024_Antarctica_Ground2', '/cresis/nvme/opr_data/accum',         'CSARP_polarimetric',        'fabric_deltak_jp'; ...
  '2023_Antarctica_Ground',  '/cresis/dataproducts/opr_data/accum', 'CSARP_polarimetric_unwrap', 'fabric_deltak'};

addpath(proj_root);                               % +ptt
addpath(fullfile(proj_root,'opr_fabric'));        % fabric_task
addpath(fullfile(proj_root,'opr_fabric','test','stubs')); % shadow opr_* helpers

t0 = tic;
n_done = 0; n_skip = 0; n_fail = 0;

for si = 1:size(product_tbl, 1)
  [season, data_root, in_name, out_path] = product_tbl{si, :};
  in_dir = fullfile(data_root, season, in_name);
  season_root = fullfile(scratch, season);
  if ~exist(season_root,'dir'), mkdir(season_root); end
  link_fn = fullfile(season_root, in_name);
  status = system(sprintf('ln -sfn ''%s'' ''%s''', in_dir, link_fn));
  if status ~= 0
    fprintf('[FAIL] %s: could not create symlink %s -> %s; skipping season\n', ...
      season, link_fn, in_dir);
    continue;
  end

  segs = dir(in_dir);
  segs = segs([segs.isdir] & ~cellfun('isempty', regexp({segs.name}, '^\d{8}_\d{2}$')));
  for gi = 1:numel(segs)
    day_seg = segs(gi).name;
    frms = dir(fullfile(in_dir, day_seg, 'Data_*.mat'));
    for fi = 1:numel(frms)
      tok = regexp(frms(fi).name, '^Data_(\d{8}_\d{2})_(\d{3})\.mat$', 'tokens');
      if isempty(tok), continue; end
      frm = str2double(tok{1}{2});

      out_fn = fullfile(season_root, ['CSARP_' out_path], day_seg, frms(fi).name);
      if exist(out_fn,'file')
        n_skip = n_skip + 1;
        continue;
      end

      param = [];
      param.day_seg = day_seg;
      param.season_name = season;
      param.radar_name = 'accum3';
      param.load.frm = frm;
      param.opr_file_lock = false;
      param.stub_out_root = season_root;

      pf = [];
      pf.in_path = in_name(7:end); % strip 'CSARP_' for the stub
      pf.out_path = out_path;
      pf.dtau_source = 'deltak';
      pf.inversion = 'joint';
      pf.reg = 0.05;
      pf.img = 0;
      pf.out_file_exts = {'.png'};
      pf.fc = 750e6;
      pf.phase_sign = -1;           % forced; settled by Ridge A cell regression (corr -0.984)
      pf.coherence_threshold = 0.5;
      pf.min_coverage = 0.3;
      pf.block_size = 1000;
      pf.num_intervals = 10;
      pf.ref_twtt_offset = 50e-9;
      pf.half_offset = 0;
      pf.ptt = struct('H', 2000, 'bco_depth', 60, 'lam_z_sfc', 1/3, 'lam_z_bed', 1/3);
      param.fabric = pf;

      t1 = tic;
      try
        ok = fabric_task(param);
        if ok
          n_done = n_done + 1;
          fprintf('[OK]   %s_%03d (%.0f s)\n', day_seg, frm, toc(t1));
        else
          n_fail = n_fail + 1;
          fprintf('[SKIP] %s_%03d (task returned false)\n', day_seg, frm);
        end
      catch ME
        n_fail = n_fail + 1;
        fprintf('[FAIL] %s_%03d: %s\n', day_seg, frm, ME.message);
      end
      close all;
    end
  end
end

fprintf('\nDelta-k batch done: %d inverted, %d already existed, %d failed/skipped (%.1f min)\n', ...
  n_done, n_skip, n_fail, toc(t0)/60);
