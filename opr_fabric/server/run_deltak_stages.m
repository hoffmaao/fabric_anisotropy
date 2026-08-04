%RUN_DELTAK_STAGES Per-stage delta-k diagnostic on a single SLC pair.
%   Runs all three dtau estimators over ONE polarimetric frame and saves
%   their depth profiles side by side, to locate where the Ridge A
%   delta-k amplitude is lost. scripts/figures/deltak_vs_joint.py shows
%   the delta-k dlam standard deviation falling from 1.3x the joint chain
%   near the surface to 18x below 1500 m; this runner exposes the
%   intermediate quantities that comparison cannot see:
%     stage A  adjacent narrow sub-band pairs (df ~ B/n_sub), smoothed
%     stage Q  quarter-band pairs, integers resolved by A
%     stage B  half-band pair, integers resolved by Q (the deliverable)
%     snaphu   ptt.blendTraveltime on the unwrapped SNAPHU phase
%     coreg    ptt.blendTraveltime on the coregistration offsets alone
%   The ladder stages share one surface-referencing convention and one
%   analysis-cell grid, so a drop in profile amplitude between A, Q and B
%   is the ladder losing signal, whereas a drop already present at A is
%   the cell averaging (dk.cell_twtt x dk.cell_ntr) low-passing dtau in
%   depth before the ladder ever runs.
%
%   Launch on mem1 with:
%     /opt/sw/matlab/2024b/bin/matlab -batch "run('<code>/opr_fabric/server/run_deltak_stages.m')"
%
%   Writes deltak_stages_<day_seg>_<frm>.mat into the scratch fabric dir,
%   for scripts/figures/deltak_stages.py.

%% Frame selection
scratch  = '/kucresis/scratch/hoffmana_sta/fabric';
season   = '2024_Antarctica_Ground2';
data_root = '/cresis/nvme/opr_data/accum';
in_name  = 'CSARP_polarimetric';
day_seg  = '20250108_02';
frm      = 1;

fc = 750e6;
phase_sign = -1;   % forced, as in run_deltak_scratch

this_dir = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(this_dir));
addpath(proj_root);                               % +ptt
addpath(fullfile(proj_root,'opr_fabric'));
addpath(fullfile(proj_root,'opr_fabric','test','stubs'));

in_fn = fullfile(data_root, season, in_name, day_seg, ...
  sprintf('Data_%s_%03d.mat', day_seg, frm));
assert(exist(in_fn,'file')==2, 'Polarimetric frame not found: %s', in_fn);

fprintf('=== delta-k stage diagnostic: %s_%03d ===\n', day_seg, frm);
fprintf('Loading %s\n', in_fn);

have = whos('-file', in_fn);
names = {have.name};
need = {'ref','sec','interferogram_coherence','Time','Surface'};
missing = setdiff(need, names);
assert(isempty(missing), 'Frame lacks %s (delta-k needs the SLCs).', strjoin(missing,', '));

want = intersect({'ref','sec','interferogram_coherence','interferogram_mlook', ...
  'snaphu_out_phase','row_offset','Time','Surface','param_polarimetric', ...
  'param_records'}, names);
pol = load(in_fn, want{:});

%% Shared map (same assembly as fabric_task)
map = [];
map.Time = pol.Time;
map.Surface = pol.Surface;
map.fc = fc;
map.coherence = abs(pol.interferogram_coherence);
% Waveform-combine boundaries, matching fabric_task's assembly so the
% seam mask behaves identically here and in the production chain
map.img_comb = [];
for pname = {'param_polarimetric','param_records'}
  p = pname{1};
  if isfield(pol,p) && isstruct(pol.(p)) && isfield(pol.(p),'array') ...
      && isfield(pol.(p).array,'img_comb') && ~isempty(pol.(p).array.img_comb)
    map.img_comb = pol.(p).array.img_comb;
    if isfield(pol.(p).array,'img_comb_mult') ...
        && ~isempty(pol.(p).array.img_comb_mult)
      map.img_comb_mult = pol.(p).array.img_comb_mult;
    end
    break;
  end
end
if isfield(pol,'row_offset')
  map.row_offset = pol.row_offset;
else
  map.row_offset = [];
end

opts = [];
opts.fc = fc;
opts.phase_sign = phase_sign;
opts.coherence_threshold = 0.5;
opts.ref_twtt_offset = 50e-9;
opts.blend_coreg_en = true;
opts.block_size = size(map.coherence,2);   % one block: whole-frame profile
opts.deltak = struct();

%% Estimator 1: delta-k (carries the per-stage cell taus in info.deltak)
map_dk = map;
map_dk.phase = [];
map_dk.phase_is_unwrapped = true;
slc = struct('ref', pol.ref, 'sec', pol.sec);
[dtau_dk, info_dk] = ptt.deltakTraveltime(slc, map_dk, opts);
clear slc;
pol = rmfield(pol, intersect({'ref','sec'}, fieldnames(pol)));
blk_dk = ptt.blockAverage(dtau_dk, map_dk, info_dk, opts);
clear dtau_dk;   % Nt x Nx doubles are ~230 MB each on an undecimated frame

%% Estimator 2: SNAPHU unwrapped phase
prof_snaphu = [];
if isfield(pol,'snaphu_out_phase')
  map_sn = map;
  map_sn.phase = pol.snaphu_out_phase;
  map_sn.phase_is_unwrapped = true;
  [dtau_sn, info_sn] = ptt.blendTraveltime(map_sn, opts);
  blk_sn = ptt.blockAverage(dtau_sn, map_sn, info_sn, opts);
  prof_snaphu = blk_sn.dtau(:,1);
  clear dtau_sn;
  fprintf('snaphu: phase_sign %+d, block fringes %d\n', ...
    info_sn.phase_sign, blk_sn.blend_fringes(1));
else
  warning('No snaphu_out_phase in this frame; the SNAPHU leg is empty.');
end

%% Estimator 3: coregistration offsets alone
prof_coreg = [];
if ~isempty(map.row_offset)
  opts_cg = opts;
  opts_cg.dtau_source = 'coreg';
  map_cg = map;
  map_cg.phase = [];
  [dtau_cg, info_cg] = ptt.blendTraveltime(map_cg, opts_cg);
  blk_cg = ptt.blockAverage(dtau_cg, map_cg, info_cg, opts_cg);
  prof_coreg = blk_cg.dtau(:,1);
  clear dtau_cg;
else
  warning('No row_offset in this frame; the coreg leg is empty.');
end

%% Reduce every rung through the identical downstream path
% Each stage is put on the full grid exactly as ptt.deltakTraveltime does
% for the dtau it ships, then block averaged with the same mask and the
% same weights. Reducing the stages one way (a cell-column statistic) and
% the phase estimators another (a coherence-weighted block average) would
% confound the reduction with the amplitude difference being measured.
stage = info_dk.deltak.stage;
t_cell = info_dk.deltak.t_cell;
x_cell = info_dk.deltak.x_cell;
Nx = size(map.coherence,2);
prof = struct();
for sname = {'tau_A','tau_Q','tau_B'}
  s = sname{1};
  g = ptt.cellToGrid(stage.(s), t_cell, x_cell, map.Time, Nx);
  b = ptt.blockAverage(g, map_dk, info_dk, opts);
  prof.(s) = b.dtau(:,1);
  clear g;
end
prof_A = prof.tau_A;
prof_Q = prof.tau_Q;
prof_B = prof.tau_B;

% Consistency check: stage B routed through cellToGrid + blockAverage here
% is the same quantity ptt.deltakTraveltime already ships as dtau, so the
% two must agree. A mismatch means the reduction in this runner has drifted
% from the estimator's own and the stage comparison below is not like-for-
% like.
d_chk = prof_B - blk_dk.dtau(:,1);
d_chk = d_chk(isfinite(d_chk));
if ~isempty(d_chk)
  fprintf('stage B vs shipped dtau: max |diff| %.3g ns (expect ~0)\n', ...
    1e9*max(abs(d_chk)));
end

Time = map.Time;
Surface = map.Surface;
prof_deltak = blk_dk.dtau(:,1);
coverage = blk_dk.coverage(:,1);
coh = blk_dk.coh(:,1);
df = info_dk.deltak.df;
margin_Q = info_dk.deltak.margin_Q;
margin_B = info_dk.deltak.margin_B;
phase_sign_dk = info_dk.phase_sign;
img_comb = map.img_comb;

%% Report the amplitude by depth band, which is what the figure quantifies
fprintf('\nProfile standard deviation by TWTT band below surface (ns):\n');
s0 = mean(Surface(isfinite(Surface)));
bands = [0 2; 2 5; 5 10; 10 20]*1e-6;
fprintf('  %-12s %8s %8s %8s %8s %8s\n','band(us)','A','Q','B','snaphu','coreg');
for bi = 1:size(bands,1)
  sel = Time(:) - s0 >= bands(bi,1) & Time(:) - s0 < bands(bi,2);
  fprintf('  %-12s %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
    sprintf('%g-%g', bands(bi,1)*1e6, bands(bi,2)*1e6), ...
    1e9*nanstd_local(prof_A(sel)), 1e9*nanstd_local(prof_Q(sel)), ...
    1e9*nanstd_local(prof_B(sel)), ...
    1e9*nanstd_local(sel_or_empty(prof_snaphu, sel)), ...
    1e9*nanstd_local(sel_or_empty(prof_coreg, sel)));
end

%% Save
out_dir = fullfile(scratch, 'stages');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
out_fn = fullfile(out_dir, sprintf('deltak_stages_%s_%03d.mat', day_seg, frm));
fprintf('\nSaving %s\n', out_fn);
save(out_fn, '-v7', 'Time', 'Surface', 't_cell', 'prof_A', 'prof_Q', 'prof_B', ...
  'prof_deltak', 'prof_snaphu', 'prof_coreg', 'coverage', 'coh', 'df', ...
  'margin_Q', 'margin_B', 'phase_sign_dk', 'img_comb', 'day_seg', 'frm');

fprintf('Done.\n');

%% ---- local helpers (Octave-safe; no stats toolbox) -------------------
function s = nanstd_local(v)
v = v(isfinite(v));
if numel(v) < 2, s = NaN; else, s = std(v); end
end

function v = sel_or_empty(p, sel)
if isempty(p), v = []; else, v = p(sel); end
end
