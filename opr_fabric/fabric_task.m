function [success] = fabric_task(param)
% [success] = fabric_task(param)
%
% Cluster task for fabric.m. Loads the CSARP_polarimetric product for one
% frame and runs the +ptt processing chain on it:
%   ptt.blendTraveltime  - traveltime differences dtau(twtt, along-track)
%                          from the interferogram phase (SNAPHU-unwrapped
%                          when available) blended with the coregistration
%                          row offsets (sign + integer fringe ambiguity)
%   ptt.blockAverage     - coherence-weighted along-track block averaging
%   ptt.invertBlocks     - inversion for the horizontal fabric contrast
%                          dlam = lam_x - lam_y through the Maxwell-Garnett
%                          firn model (Rathmann 2026), by exact layer
%                          stripping or the smoothness-regularized joint
%                          solve selected by param.fabric.inversion
% All processing options are taken from the param.fabric struct (see
% fabric.m), which the +ptt functions read directly.
%
% Sign/axis convention: the polarimetric product's interferogram is
% sec .* conj(ref) (ref = HH and sec = VV rotated by synth_rot_deg), and
% coregistration row_offset > 0 means sec arrives later. dtau here is
% t_sec - t_ref, so dlam > 0 means the c-axis concentration is greater
% along the secondary (rotated V) axis than the reference (rotated H) axis.
%
% param = struct controlling the loading and fabric processing

% er_ice, c = speed of light
physical_constants;
[output_dir,radar_type,radar_name] = opr_output_dir(param.radar_name);

fc = param.fabric.fc;

frm = param.load.frm;
frm_id = sprintf('%s_%03d', param.day_seg, frm);
fprintf('fabric processing frame %s (%s)\n', frm_id, datestr(now));

% fn_name: data filename string Data_YYYYMMDD_SS_FFF.mat or
%             Data_img_II_YYYYMMDD_SS_FFF.mat
if param.fabric.img == 0
  fn_name = sprintf('Data_%s.mat', frm_id);
else
  fn_name = sprintf('Data_img_%02d_%s.mat', param.fabric.img, frm_id);
end

out_dir = opr_filename_out(param,param.fabric.out_path);
if ~exist(out_dir,'dir')
  mkdir(out_dir);
end
out_fn = fullfile(out_dir,fn_name);
[~,fn_name] = fileparts(fn_name);

in_fn_dir = opr_filename_out(param,param.fabric.in_path);
in_fn = fullfile(in_fn_dir,[fn_name '.mat']);
if ~exist(in_fn,'file')
  warning('The polarimetric file does not exist for this frame. Skipping this frame. Perhaps param.fabric.in_path is incorrect or polarimetric.m has not been run. File does not exist:\n  %s.', in_fn);
  success = false;
  return;
end

%% Load polarimetric product and assemble the input map
% =========================================================================
% Load only the fields the chain uses. The product also carries the full
% complex ref/sec/sec_reg images (GBs in standardphase products): the
% phase/coreg modes never touch them, but dtau_source = 'deltak' needs
% the ref and UNREGISTERED sec SLCs (coregistration shifts the envelope,
% which removes exactly the group delay delta-k measures).
deltak_en = isfield(param.fabric,'dtau_source') ...
  && strcmp(param.fabric.dtau_source,'deltak');
if deltak_en
  want = {'interferogram_coherence','row_offset','Time','GPS_time', ...
    'Latitude','Longitude','Elevation','Surface','param_records', ...
    'param_polarimetric','ref','sec'};
  required = {'Time','Surface','interferogram_coherence','ref','sec'};
else
  want = {'interferogram_mlook','interferogram_coherence','snaphu_out_phase', ...
    'row_offset','Time','GPS_time','Latitude','Longitude','Elevation', ...
    'Surface','param_records','param_polarimetric'};
  required = {'Time','Surface','interferogram_coherence','interferogram_mlook'};
end
have = whos('-file', in_fn);
missing = setdiff(required, {have.name});
if ~isempty(missing)
  warning('Required variable(s) %s missing from polarimetric file. Skipping this frame. Perhaps param.fabric.in_path points at the wrong product. File:\n  %s.', strjoin(missing, ', '), in_fn);
  success = false;
  return;
end
sel = intersect(want, {have.name});
pol = load(in_fn, sel{:});

map = [];
map.Time = pol.Time;
map.Surface = pol.Surface;
map.fc = fc;
map.coherence = abs(pol.interferogram_coherence);

% Waveform-combine boundaries, so ptt.surfaceReference can drop the seams
% (ptt.imgCombSeam). Prefer the polarimetric product's own array settings:
% they describe how the images this interferogram was formed from were
% stitched. param_records is the fallback for products that did not carry
% their array params forward. img_comb_mult, the surface multiplier
% img_combine.m applies on top of img_comb, rides along when the product
% carries it (ptt.imgCombSeam defaults it to Inf otherwise, as OPR does).
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

if deltak_en
  % Delta-k derives dtau from the SLC spectra directly; no phase field
  map.phase = [];
  map.phase_is_unwrapped = true;
elseif param.fabric.use_snaphu_phase && isfield(pol,'snaphu_out_phase')
  map.phase = pol.snaphu_out_phase;
  map.phase_is_unwrapped = true;
else
  if param.fabric.use_snaphu_phase
    warning('snaphu_out_phase not found in polarimetric product; falling back to wrapped interferogram phase. Fringe ambiguities will be resolved per pixel from the coregistration row offsets, which limits accuracy to the coregistration noise.');
  end
  map.phase = angle(pol.interferogram_mlook);
  map.phase_is_unwrapped = false;
end

if isfield(pol,'row_offset')
  % Coregistration row_offset > 0: sec pixel belongs row_offset rows above
  % its current location, i.e. sec arrives later by row_offset*dt
  map.row_offset = pol.row_offset;
else
  map.row_offset = [];
  warning('row_offset not found in polarimetric product (coregistration disabled?). Coregistration blending is unavailable; the phase sign cannot be auto-detected.');
end

%% Traveltime differences, block averaging, inversion
% =========================================================================
if deltak_en
  slc = struct('ref', pol.ref, 'sec', pol.sec);
  % Drop the product's own copies: pol stays live through the inversion,
  % plotting and save, and the SLCs are GB-scale on undecimated frames
  pol = rmfield(pol, {'ref','sec'});
  [dtau, info] = ptt.deltakTraveltime(slc, map, param.fabric);
  clear slc;
else
  [dtau, info] = ptt.blendTraveltime(map, param.fabric);
end
phase_sign = info.phase_sign;
fprintf('Phase sign: %+d\n', phase_sign);

blk = ptt.blockAverage(dtau, map, info, param.fabric);
Nblk = numel(blk.starts);

GPS_time = cellfun(@(c) mean(pol.GPS_time(c)), blk.cols);
Latitude = cellfun(@(c) mean(pol.Latitude(c)), blk.cols);
Longitude = cellfun(@(c) mean(pol.Longitude(c)), blk.cols);
Elevation = cellfun(@(c) mean(pol.Elevation(c)), blk.cols);
Surface = blk.Surface;

% Firn-ice column model
par = ptt.defaultParams();
for fld = fieldnames(param.fabric.ptt).'
  par.(fld{1}) = param.fabric.ptt.(fld{1});
end
if isfield(par,'bco_depth')
  if ~isempty(par.bco_depth)
    par.zhat_bco = 1 - par.bco_depth/par.H;
  end
  par = rmfield(par,'bco_depth');
end

inv = ptt.invertBlocks(blk, map, par, param.fabric);

dlam = inv.dlam;
dlam_top_depth = inv.top_depth;
dlam_bot_depth = inv.bot_depth;
dtau_obs = inv.dtau_obs;
dtau_fit = inv.dtau_fit;
dlam_quality = inv.quality;
% 1 where the interval's dtau observation was interpolated across a gap
% (waveform-combine seam, incoherent run) rather than measured there;
% dlam_quality cannot show this, it is the unmasked coherence
dlam_interpolated = inv.interpolated;
dlam_clipped = inv.clipped;
% 1 where the interval's dlam shares its information with the reference
% offset (the shallowest interval when obs.zref is in play): quotable
% fabric starts below it. ref_offset is the recovered reference error [ns].
dlam_ref_degenerate = inv.ref_degenerate;
ref_offset = inv.ref_offset;
dtau_rms = inv.rms;
reg_alpha = inv.alpha;
dtau_blk = blk.dtau;
coh_blk = blk.coh;
coverage_blk = blk.coverage;
blend_fringes = blk.blend_fringes;

%% Plot results
% =========================================================================
h_fig = get_figures(2,true);
clear h_axes;
if exist('parula','file') || exist('parula','builtin')
  cmap = parula(256);
else
  cmap = jet(256); % Octave fallback
end

% Docking is unavailable in headless batch sessions (-nodisplay)
dock_en = usejava('desktop');
fig_idx = 1; clf(h_fig(fig_idx));
if dock_en, set(h_fig(fig_idx),'WindowStyle','docked'); end
h_axes(fig_idx) = axes('parent',h_fig(fig_idx));
mean_bot_depth = mean(dlam_bot_depth,2,'omitnan');
imagesc(1:Nblk, mean_bot_depth, dlam, 'parent', h_axes(fig_idx));
set(h_axes(fig_idx),'YDir','reverse');
colormap(h_axes(fig_idx),cmap);
h_colorbar = colorbar(h_axes(fig_idx));
set(get(h_colorbar,'ylabel'),'string','\Delta\lambda = \lambda_{sec} - \lambda_{ref}')
title(h_axes(fig_idx),sprintf('Horizontal fabric contrast %s',regexprep(frm_id,'_','\\_')));
xlabel(h_axes(fig_idx),'Block');
ylabel(h_axes(fig_idx),'Interval bottom depth (m)');
grid(h_axes(fig_idx),'on');

fig_idx = 2; clf(h_fig(fig_idx));
if dock_en, set(h_fig(fig_idx),'WindowStyle','docked'); end
h_axes(fig_idx) = axes('parent',h_fig(fig_idx));
% Draw only the cells the inversion uses: below min_coverage the block
% average is a handful of pixels snapped to whole fringes, which drew as
% horizontal jumps at the bottom of the record
dtau_plot = dtau_blk;
dtau_plot(coverage_blk < param.fabric.min_coverage) = NaN;
plot(h_axes(fig_idx), 1e9*dtau_plot, (pol.Time - mean(Surface,'omitnan'))*1e6);
set(h_axes(fig_idx),'YDir','reverse');
title(h_axes(fig_idx),sprintf('Block-averaged \\Delta\\tau %s',regexprep(frm_id,'_','\\_')));
xlabel(h_axes(fig_idx),'\Delta\tau (ns)');
ylabel(h_axes(fig_idx),'TWTT below surface (\mus)');
grid(h_axes(fig_idx),'on');

%% Save outputs
% =========================================================================
Time = pol.Time;
if isfield(pol,'param_records')
  param_records = pol.param_records;
else
  param_records = [];
end
if isfield(pol,'param_polarimetric')
  param_polarimetric = pol.param_polarimetric;
else
  param_polarimetric = [];
end
param_fabric = param;
if param.opr_file_lock
  file_version = '1L';
else
  file_version = '1';
end
file_type = 'fabric';

fprintf('Saving output file:\n  %s\n', out_fn);
opr_save(out_fn,'dlam','dlam_top_depth','dlam_bot_depth','dlam_quality', ...
  'dlam_interpolated','dlam_clipped','dlam_ref_degenerate','ref_offset', ...
  'dtau_rms','reg_alpha', ...
  'dtau_obs','dtau_fit','dtau_blk','coh_blk','coverage_blk','blend_fringes', ...
  'phase_sign','Time','GPS_time','Latitude','Longitude','Elevation','Surface', ...
  'param_fabric','param_polarimetric','param_records','file_type','file_version');

for file_ext = param.fabric.out_file_exts
  file_ext = file_ext{1};
  fprintf('Saving output images of type %s\n', file_ext);
  opr_saveas(h_fig(1),fullfile(out_dir,sprintf('%s_dlam%s',fn_name,file_ext)));
  opr_saveas(h_fig(2),fullfile(out_dir,sprintf('%s_dtau%s',fn_name,file_ext)));
end

%% Done
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
