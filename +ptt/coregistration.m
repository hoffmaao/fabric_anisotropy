function [img2_coregistered,row_offset_full,col_offset_full,tiles] = coregistration(img1,img2,options)
% [img2_coregistered,row_offset_full,col_offset_full,tiles] = ptt.coregistration(img1,img2,options)
%
% Performs 2D coregistration
%
% PROVENANCE. A copy of the CReSIS/OPR toolbox function
% opr/matlab/processing/coregistration.m (toolbox commit ddb056a, 10 Sep
% 2025), brought into this package on 2 Sep 2026 so it can be MODIFIED for
% polarimetric channel pairs. ptt.coregisterChannels deliberately called
% the toolbox version until now, so the quad-pol path stayed aligned by the
% same code as the shipped CSARP_polarimetric products; that reasoning
% still holds wherever this copy is run with its toolbox defaults, and
% every change below is switchable so the toolbox behaviour is one option
% away. The toolbox function keeps its name; this one lives in +ptt and
% never shadows it.
%
% WHY A MODIFIED COPY, AND WHAT WAS FOUND. The copy was made to test a
% hypothesis about the Thwaites margin, where HH-VV coherence dies and the
% co-pol phase stops changing with depth: that each channel there is a
% mixture of both birefringent eigenmodes, so the tile cross-correlation
% has peaks at lag 0 and at +-dtau and an independent per-tile argmax
% flips between them. MEASURED on the shipped product of 20240108_01_001
% (scripts/prototypes/coreg_lag_diag.py, 2 Sep 2026): NOT SO. The HH-VV
% correlation has ONE peak, drifting only 1-3 bins over 1500 m, and where
% coherence is dead it is dead at every lag within +-12 bins (+-40 ns),
% so no shift recovers it; registration moves band-median |C| by at most
% +0.09 and makes 4.8% of cells worse. The loss is at the single-look
% level (true |C| below ~0.2 at 40-60 dB SNR in the western fringe fan:
% the pair is depolarised, not misregistered) and, below ~1400 m, SNR.
% What this copy DOES fix is a real defect that the same product shows:
% offsets of -21.8 bins from a +-5 search, the vertex of a parabola
% fitted to a flat correlation, applied verbatim by the resampler
% (70% of cells below 1500 m). See MODIFICATIONS in the options block;
% both are off by default so the toolbox behaviour is one option away.
%
% img1 = image 1 (master or reference image) matrix of size Nt by Nx
%
% img2 = image 2 (slave or secondary image) matrix of size Nt by Nx
%
% img2_coregistered = img2 sinc interpolated to be registered to img1. Matrix matches img1 size.
%
% row_offset_full = number of rows that img2 needs to be shifted by. Matrix
% matches img1 size. For example, a value of 2 means that the corresponding
% pixel in img2 belongs in a pixel 2 rows above its current location.
% 
% col_offset_full = same as row_offset_full except for column offsets
%
% tiles: structure with debugging fields
arguments
  img1 {mustBeNumeric}
  img2 {mustBeNumeric}
  options.Tt (1,1) {mustBeNumeric} = 31
  options.Tx (1,1) {mustBeNumeric} = 31
  options.overlap_t (1,1) {mustBeNumeric} = 24
  options.overlap_x (1,1) {mustBeNumeric} = 24
  options.search_t (1,1) {mustBeNumeric} = 9
  options.search_x (1,1) {mustBeNumeric} = 9
  options.Mt (1,1) {mustBeNumeric} = 4
  options.Mx (1,1) {mustBeNumeric} = 4
  options.xcorr_t (1,1) {mustBeNumeric} = 2
  options.xcorr_x (1,1) {mustBeNumeric} = 2
  options.row_offset_medfilt_t (1,1) {mustBeNumeric} = 1 % +/- this many samples used in medfilt1
  options.row_offset_medfilt_x (1,1) {mustBeNumeric} = 1
  options.col_offset_medfilt_t (1,1) {mustBeNumeric} = 1
  options.col_offset_medfilt_x (1,1) {mustBeNumeric} = 1
  options.xcorr_trim_t (1,1) {mustBeNumeric} = 1 % This many samples trimmed after medfilt1
  options.xcorr_trim_x (1,1) {mustBeNumeric} = 1
  options.power_detect_en (1,1) {mustBeNumericOrLogical} = false
  options.one_dim_search_en (1,1) {mustBeNumericOrLogical} = false
  options.debug_plot_en (1,1) {mustBeNumericOrLogical} = false
  options.debug_print_en (1,1) {mustBeNumericOrLogical} = true
  % --- MODIFICATIONS (all off by default = toolbox behaviour) ---
  % peak_clamp_en: keep the quadratic sub-sample refinement inside the
  % +-xcorr_t patch it was fitted on. The vertex of a parabola fitted to a
  % flat or noisy correlation lands anywhere (a measured -21.8 bins on a
  % Thwaites margin product whose search was +-5), and the sinc resampler
  % then applies that shift verbatim. A refinement that leaves its patch is
  % not a refinement; the integer peak is kept instead.
  options.peak_clamp_en (1,1) {mustBeNumericOrLogical} = false
  % peak_contrast_min: tiles whose correlation peak is not at least this
  % factor (linear) above the median of the correlation over the search
  % range carry no offset information; their offset is marked NaN and
  % filled from neighbours by the interpolation instead of voting.
  options.peak_contrast_min (1,1) {mustBeNumeric} = 0
end

%% Input and Setup

% Tile size, must be odd integer, measured in pixels
Tt = options.Tt;
Tx = options.Tx;
% Overlap tile, must be positive integer, measured in pixels
overlap_t = options.overlap_t;
overlap_x = options.overlap_x;
% Search range, must be nonegative integer, measured in pixels
search_t = options.search_t;
search_x = options.search_x;
% Oversampling rate for cross correlation
Mt = options.Mt;
Mx = options.Mx;
% Quadratic interpolation patch size (+/- this many samples used)
xcorr_t = options.xcorr_t;
xcorr_x = options.xcorr_x;
% Enable/disable power detection before cross correlation
power_detect_en = options.power_detect_en;
% One dimensional search for correlation offset (i.e. col_offset_full will
% be zero)
one_dim_search_en = options.one_dim_search_en;

Nt = size(img1,1);
Nx = size(img1,2);

img1(~isfinite(img1)) = 0;

img2(~isfinite(img2)) = 0;

Nr_M = Nt*Mt;
Nc_M = Nx*Mx;
f_row_M = 1/Mt * (-floor(Nr_M/2) : floor((Nr_M-1)/2)).';
f_col_M = 1/Mx * -floor(Nc_M/2) : floor((Nc_M-1)/2);

%% xcorr of each tile
cols = 1:Tx-overlap_x:Nx-Tx+1;
rows = 1:Tt-overlap_t:Nt-Tt+1;
offset_t = zeros(length(rows),length(cols));
offset_x = zeros(length(rows),length(cols));
peak_val = zeros(length(rows),length(cols));
total_val = zeros(length(rows),length(cols));
contrast = inf(length(rows),length(cols));

start_time = tic;
debug_plot_en = options.debug_plot_en;
for col_idx = 1:length(cols)
  % Estimate the amount of time left to complete this loop
  time_left_estimate_sec = toc(start_time)/(col_idx-1)*(length(cols)-col_idx+1);
  if options.debug_print_en
    fprintf('%s: %d of %d: %.0f sec elapsed, %.0f sec left\n', mfilename, col_idx, length(cols), toc(start_time), time_left_estimate_sec);
  end

  col = cols(col_idx);

  for row_idx = 1:length(rows)
    row = rows(row_idx);

    % HACK TO DEBUG AND TEST
    % img2 = img1;
    % row = round(rows(end)/2); col = round(cols(end)/2);

    % Create tiles
    % =====================================================================
    tile1 = zeros(2*(Tt+2*search_t)-1,2*(Tx+2*search_x)-1);
    tile2 = zeros(2*(Tt+2*search_t)-1,2*(Tx+2*search_x)-1);

    tile1(1:Tt,1:Tx) = img1(row:row+Tt-1,col:col+Tx-1);

    row_first = max(1,search_t-row+1);
    row_last = Tt+search_t*2-1-max(0,row+Tt+search_t-Nt);

    col_first = max(1,search_x-col+1);
    col_last = Tx+search_x*2-1-max(0,col+Tx+search_x-Nx);

    tile2(row_first:row_last,col_first:col_last) ...
      = img2(row-search_t+(row_first:row_last),col-search_x+(col_first:col_last));

    if debug_plot_en
      % Debug: Plot tiles so they are registered with each other (i.e. they
      % should line up perfectly if the images are already co-registered)
      figure(1); clf;
      imagesc(search_t-1 + (1:size(tile1,2)), search_x-1 + (1:size(tile1,1)), db(tile1));
      figure(2); clf;
      imagesc(db(tile2));
      link_figures([1 2]);
    end

    % Magnitude only correlation
    % =====================================================================
    if power_detect_en
      tile1 = abs(tile1).^2;
      tile2 = abs(tile2).^2;
    end

    % Cross correlation in the frequency domain
    % =====================================================================
    xtile_raw = ifft2(conj(fft2(tile1-mean(tile1(:)))) .* fft2(tile2-mean(tile2(:))));

    % Overinterpolate cross correlation result in the frequency domain
    % =====================================================================
    xtile = interpft(interpft(xtile_raw,Mt*size(xtile_raw,1),1),Mx*size(xtile_raw,2),2);
    xtile = fftshift(xtile);

    % Find the peak
    % =====================================================================
    if one_dim_search_en
      [max_vals,max_rows] = max(xtile);
      max_col = floor(size(xtile,2)/2) + 1 + Mx*(search_x-1);
      max_val = max_vals(max_col);
      max_row = max_rows(max_col);

    else
      % Two-dimensional search
      [max_vals,max_rows] = max(xtile);
      [max_val,max_col] = max(max_vals);
      max_row = max_rows(max_col);
    end

    if debug_plot_en
      % Debug: Plot tiles
      figure(3); clf;
      imagesc(db(xtile_raw));
      figure(4); clf;
      imagesc(db(xtile));
      hold on
      plot(max_col,max_row,'kx');
    end

    if one_dim_search_en
      % 1D Quadratic interpolation centered on the maximum correlation
      % 1 + y + y^2
      % =====================================================================
      A = [];
      b = [];
      % Test
      %C = randn(6,1)
      for row_rel = -xcorr_t:xcorr_t
        A(end+1,:) = [1 row_rel row_rel^2];
        b(end+1,1) = db(xtile(max_row+row_rel,max_col));
        % Test
        %b(end,1) = A(end,:)*C;
      end
      % Least squares solution to the 2D quadratic fit to our
      % cross-correlation peak
      C = pinv(A)*b;

      % Solution for the stationary point of the 2D quadratic function allows
      % us to find the row and column of the peak value of the quadratic
      max_row_rel = -C(2)/(2*C(3));
      if options.peak_clamp_en && (~isfinite(max_row_rel) || C(3) >= 0 ...
          || abs(max_row_rel) > xcorr_t)
        % no interior maximum inside the patch: keep the integer peak
        max_row_rel = 0;
      end

      % Store the peak_val for debugging and evaluation purposes
      peak_val(row_idx,col_idx) = [1 max_row_rel max_row_rel^2]*C;
      total_val(row_idx,col_idx) = db(sum(abs(xtile).^2,'all'),'power');
      if options.peak_contrast_min > 0
        % contrast of the peak against the search range it was picked from
        rr = max(1, max_row - Mt*search_t) : min(size(xtile,1), max_row + Mt*search_t);
        contrast(row_idx,col_idx) = abs(xtile(max_row,max_col)) / max(median(abs(xtile(rr,max_col))), realmin);
      end

      % Convert relative values to absolute values
%          offset_x(row_idx,col_idx) = (max_col_int - floor(size(xtile,2)/2) - 1)/Mx - search_x+1;

      max_col_int = max_col;
      max_row_int = max_row + max_row_rel;

    else
      % 2D Quadratic interpolation centered on the maximum correlation
      % 1 + x + y + x*y + x^2 + y^2
      % =====================================================================
      A = [];
      b = [];
      % Test
      %C = randn(6,1)
      for col_rel = -xcorr_x:xcorr_x
        for row_rel = -xcorr_t:xcorr_t
          A(end+1,:) = [1 col_rel row_rel col_rel*row_rel col_rel^2 row_rel^2];
          b(end+1,1) = db(xtile(max_row+row_rel,max_col+col_rel));
          % Test
          %b(end,1) = A(end,:)*C;
        end
      end
      % Least squares solution to the 2D quadratic fit to our
      % cross-correlation peak
      C = pinv(A)*b;

      % Solution for the stationary point of the 2D quadratic function allows
      % us to find the row and column of the peak value of the quadratic
      max_row_rel = (2*C(5)*C(3) - C(2)*C(4)) / (C(4)^2 - 4*C(5)*C(6));
      max_col_rel = (-C(2) - C(4)*max_row_rel) / (2*C(5));
      if options.peak_clamp_en && (~all(isfinite([max_row_rel max_col_rel])) ...
          || abs(max_row_rel) > xcorr_t || abs(max_col_rel) > xcorr_x)
        max_row_rel = 0; max_col_rel = 0;
      end

      % Store the peak_val for debuggin and evaluation purposes
      peak_val(row_idx,col_idx) = [1 max_col_rel max_row_rel max_col_rel*max_row_rel max_col_rel^2 max_row_rel^2]*C;
      total_val(row_idx,col_idx) = sum(abs(xtile).^2,'all');

      % Convert relative values to absolute values
      max_col_int = max_col + max_col_rel;
      max_row_int = max_row + max_row_rel;
    end

    % Save offsets
    % =====================================================================
    offset_t(row_idx,col_idx) = (max_row_int - floor(size(xtile,1)/2) - 1)/Mt - search_t+1;
    offset_x(row_idx,col_idx) = (max_col_int - floor(size(xtile,2)/2) - 1)/Mx - search_x+1;

    if debug_plot_en
      %% xcorr of each tile - Debug Plots
      figure(5); clf;
      %       imagesc(db(xtile) - db(max_val));
      %       caxis([-6 0]);
      imagesc(-xcorr_x:xcorr_x, -xcorr_t:xcorr_t, db(xtile(max_row+(-xcorr_x:xcorr_x),max_col+(-xcorr_t:xcorr_t))))
      colorbar
      hold on;
      plot(max_col-max_col, max_row-max_row,'bx','linewidth',2);
      plot(max_col_int-max_col, max_row_int-max_row,'kx','linewidth',2);

      figure(6); clf;
      %       imagesc(db(xtile) - db(max_val));
      %       caxis([-6 0]);
      if one_dim_search_en
        imagesc(0, -xcorr_t:xcorr_t, ...
          reshape(A*C, length(-xcorr_t:xcorr_t), 1) )
      else
        imagesc(-xcorr_x:xcorr_x, -xcorr_t:xcorr_t, ...
          reshape(A*C, length(-xcorr_t:xcorr_t), length(-xcorr_x:xcorr_x)) )
      end
      colorbar
      hold on;
      plot(max_col-max_col, max_row-max_row,'bx','linewidth',2);
      plot(max_col_int-max_col, max_row_int-max_row,'kx','linewidth',2);
      title(sprintf('(%d,%d) is (%.2f,%.2f)', row, col, offset_t(row_idx,col_idx), offset_x(row_idx,col_idx)));
      link_figures([5 6]);
    end

    if 0
      % Debug to print out offsets
      offset_t(row_idx,col_idx)
      offset_x(row_idx,col_idx)
    end
  end

end

%% Smooth xcorr offsets
% =========================================================================
if options.peak_contrast_min > 0
  % a tile with no usable peak does not vote: its offset is replaced by the
  % median of its finite neighbours before the median filter runs, so it
  % can neither pull the filter nor reach the interpolant as a value
  bad = contrast < options.peak_contrast_min;
  offset_t(bad) = NaN;
  offset_x(bad) = NaN;
  offset_t = H_fill_nan(offset_t);
  offset_x = H_fill_nan(offset_x);
  if options.debug_print_en
    fprintf('%s: %d of %d tiles below peak contrast %.2f, filled from neighbours\n', ...
      mfilename, nnz(bad), numel(bad), options.peak_contrast_min);
  end
end
row_offsets = medfilt2(offset_t, ...
  [2*options.row_offset_medfilt_t-1, 2*options.row_offset_medfilt_x-1],'symmetric');
col_offsets = medfilt2(offset_x, ...
  [2*options.col_offset_medfilt_t-1, 2*options.col_offset_medfilt_x-1],'symmetric');
row_offsets = row_offsets(1+options.xcorr_trim_t:end-options.xcorr_trim_t, ...
  1+options.xcorr_trim_x:end-options.xcorr_trim_x);
col_offsets = col_offsets(1+options.xcorr_trim_t:end-options.xcorr_trim_t, ...
  1+options.xcorr_trim_x:end-options.xcorr_trim_x);

if 0
  % Plot results
  figure(12); clf;
  imagesc(rows,cols,row_offsets);
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Offsets (rows)');
  title('Row offsets');
  xlabel('Row');
  ylabel('Column');

  figure(13); clf;
  imagesc(rows,cols,col_offsets);
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Offsets (columns)');
  title('Column offsets');
  xlabel('Row');
  ylabel('Column');

  link_figures([12 13]);
end

tiles = struct();
tiles.rows = rows;
tiles.cols = cols;
tiles.offset_t = offset_t;
tiles.offset_x = offset_x;
tiles.peak_val = peak_val;
tiles.total_val = total_val;
tiles.contrast = contrast;

%% Interpolate sparse xcorr results to every pixel
% =========================================================================

% Create interpolation function using scatteredInterpolant.m
[row_mesh,col_mesh] = meshgrid(cols(1+options.xcorr_trim_x:end-options.xcorr_trim_x), ...
  rows(1+options.xcorr_trim_t:end-options.xcorr_trim_t));
Frow = scatteredInterpolant(row_mesh(:)+Tt/2+1,col_mesh(:)+Tx/2+1,row_offsets(:),'natural','nearest');
Fcol = scatteredInterpolant(row_mesh(:)+Tt/2+1,col_mesh(:)+Tx/2+1,col_offsets(:),'natural','nearest');

% Interpolate
[row_mesh,col_mesh] = meshgrid(1:size(img1,2),1:size(img1,1));
row_offset_full = Frow(row_mesh,col_mesh);
col_offset_full = Fcol(row_mesh,col_mesh);

if options.debug_print_en
  fprintf('%s: Offsets interpolation complete (elapsed time %.0f sec)\n', mfilename, toc(start_time));
end

if 0
  % Plot results
  figure(12); clf;
  imagesc(rows,cols,row_offset_full);
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Offsets (rows)');
  title('Row offsets interpolated');
  xlabel('Row');
  ylabel('Column');

  figure(13); clf;
  imagesc(rows,cols,col_offset_full);
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Offsets (columns)');
  title('Column offsets interpolated');
  xlabel('Row');
  ylabel('Column');
end

%% Re-sample Image
% =========================================================================

% Preallocate registered image
img2_coregistered = zeros(size(img1));

% Create sinc interpolation window weights
Nr_sinc = 4;
Nc_sinc = 4;
Hwind = kaiser(1+2*Nr_sinc)*kaiser(1+2*Nc_sinc).';

% Iterate through every pixel in the resampled image
cols_all = 1:size(img1,2);
rows_all = 1:size(img1,1);
rel_rows = (-Nr_sinc:Nr_sinc).';
rel_cols = -Nc_sinc:Nc_sinc;
for col_idx = 1:length(cols_all)
  col = cols_all(col_idx);
  for row_idx = 1:length(rows_all)
    row = rows_all(row_idx);

    % Grab the coregistration offset for this pixel
    row_off = -row_offset_full(row_idx,col_idx);
    col_off = -col_offset_full(row_idx,col_idx);

    % Convert offset into integer offset and fractional offset
    row_center = round(row_off);
    row_frac = row_off - row_center;
    col_center = round(col_off);
    col_frac = col_off - col_center;

    rel_row_idxs = max(1,1+Nr_sinc-(row-row_center)+1) : min(2*Nr_sinc+1,1+Nr_sinc+(length(rows_all)-(row-row_center)));
    rel_col_idxs = max(1,1+Nc_sinc-(col-col_center)+1) : min(2*Nc_sinc+1,1+Nc_sinc+(length(cols_all)-(col-col_center)));
    % Apply 2D sinc interpolation centered on the target image's pixel offset
    if row-row_center >= 1 && row-row_center <= length(rows_all) ...
        && col-col_center >= 1 && col-col_center <= length(cols_all)
      img2_coregistered(row_idx,col_idx) = sum(sum(img2(row-row_center+rel_rows(rel_row_idxs),col-col_center+rel_cols(rel_col_idxs)) ...
        .* Hwind(1+Nr_sinc+rel_rows(rel_row_idxs),1+Nc_sinc+rel_cols(rel_col_idxs)) ...
        .* bsxfun(@times,sinc(row_frac+rel_rows(rel_row_idxs)),sinc(col_frac+rel_cols(rel_col_idxs))) ));
    end

  end
end
if options.debug_print_en
  fprintf('%s: Coregistration complete (elapsed time %.0f sec)\n', mfilename, toc(start_time));
end

if 0
  %% Debug plot results
  figure(14); clf;
  imagesc(db(img1));
  h_coreg_axes = gca;
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Relative power (dB)');
  title('Reference SAR Image');
  xlabel('Column');
  ylabel('Row');

  figure(15); clf;
  imagesc(db(img2));
  h_coreg_axes(2) = gca;
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Relative power (dB)');
  title('Target SAR Image');
  xlabel('Column');
  ylabel('Row');

  figure(16); clf;
  imagesc(db(img2_coregistered));
  h_coreg_axes(3) = gca;
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Relative power (dB)');
  title('Target SAR Image');
  xlabel('Column');
  ylabel('Row');

  linkaxes(h_coreg_axes)
end

function A = H_fill_nan(A)
% replace NaN entries by the median of their finite 3x3 neighbours,
% repeating until none remain (or nothing finite is left to fill from)
for it = 1:numel(A)
  bad = ~isfinite(A);
  if ~any(bad(:)), return; end
  Ap = padarray(A, [1 1], NaN);
  filled = false;
  [ib, jb] = find(bad);
  for k = 1:numel(ib)
    nb = Ap(ib(k):ib(k)+2, jb(k):jb(k)+2);
    nb = nb(isfinite(nb));
    if ~isempty(nb)
      A(ib(k), jb(k)) = median(nb);
      filled = true;
    end
  end
  if ~filled, return; end
end
