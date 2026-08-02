%RUN_FABRIC_SCRATCH Batch fabric inversion on CReSIS-server products.
%   Runs fabric_task over every frame of the CSARP_polarimetric_unwrap
%   products (Hoffman/Christianson runs, 2022 + 2023 Antarctica Ground
%   seasons), writing CSARP_fabric outputs to the user's scratch instead of
%   the shared season tree. Inputs are reached through per-season symlinks
%   created under the scratch root, so the opr_* path stubs (one shared
%   root for in/out) work unchanged.
%
%   Rerun-safe: frames whose output file already exists are skipped.
%
%   Launch on mem1 with:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('<code>/opr_fabric/server/run_fabric_scratch.m')"

data_root = '/cresis/dataproducts/opr_data/accum';
scratch   = '/kucresis/scratch/hoffmana_sta/fabric';
seasons   = {'2022_Antarctica_Ground', '2023_Antarctica_Ground'};
in_name   = 'CSARP_polarimetric_unwrap';

this_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));
addpath(proj_root);                               % +ptt
addpath(fullfile(proj_root,'opr_fabric'));        % fabric_task
addpath(fullfile(proj_root,'opr_fabric','test','stubs')); % shadow opr_* helpers

t0 = tic;
n_done = 0; n_skip = 0; n_fail = 0;

for si = 1:numel(seasons)
  season = seasons{si};
  in_dir = fullfile(data_root, season, in_name);
  season_root = fullfile(scratch, season);
  if ~exist(season_root,'dir'), mkdir(season_root); end
  % Symlink the input product into the scratch season root
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
      if isempty(tok), continue; end % skip figure sidecars etc.
      frm = str2double(tok{1}{2});

      out_fn = fullfile(season_root, 'CSARP_fabric', day_seg, frms(fi).name);
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
      pf.out_path = 'fabric';
      pf.img = 0;
      pf.out_file_exts = {'.png'};
      pf.fc = 750e6;
      pf.dtau_source = 'phase';
      pf.use_snaphu_phase = true;
      pf.blend_coreg_en = true;
      pf.phase_sign = 0;
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

fprintf('\nBatch done: %d inverted, %d already existed, %d failed/skipped (%.1f min)\n', ...
  n_done, n_skip, n_fail, toc(t0)/60);
