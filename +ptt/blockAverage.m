function blk = blockAverage(dtau, map, info, opts)
%BLOCKAVERAGE Coherence-weighted along-track block averaging of dtau.
%   blk = BLOCKAVERAGE(dtau, map, info, opts) averages the traveltime
%   difference map (from ptt.blendTraveltime) into along-track blocks of
%   opts.block_size columns, weighting by coherence^2 within the coherence
%   mask. For unwrapped phase with coregistration available, one integer
%   fringe constant per block is estimated from the median disagreement
%   between the two estimators and applied (corrects unwrapping errors and
%   any residual constant).
%
%   map, info: as produced for/by ptt.blendTraveltime.
%   opts fields (optional): block_size (1000), blend_coreg_en (true).
%
%   blk fields (Nblk blocks):
%     starts (1 x Nblk), cols {1 x Nblk} column index ranges
%     dtau (Nt x Nblk) [s], coverage (Nt x Nblk) coherent fraction,
%     coh (Nt x Nblk) mean coherence, Surface (1 x Nblk) [s],
%     blend_fringes (1 x Nblk)

if ~isfield(opts,'block_size') || isempty(opts.block_size)
  opts.block_size = 1000;
end
if ~isfield(opts,'blend_coreg_en') || isempty(opts.blend_coreg_en)
  opts.blend_coreg_en = true;
end

Nt = size(dtau,1);
Nx = size(dtau,2);

% NaN coherence or dtau samples must carry zero weight AND zero value:
% NaN*0 = NaN in MATLAB, so either would otherwise NaN the whole row of
% every block they fall in despite the mask.
w = map.coherence.^2 .* info.coh_mask;
bad = ~isfinite(w) | ~isfinite(dtau);
w(bad) = 0;
dtau_w = dtau;
dtau_w(bad) = 0;

blk.starts = 1:opts.block_size:Nx;
Nblk = numel(blk.starts);
blk.cols = cell(1,Nblk);
blk.dtau = nan(Nt,Nblk);
blk.coverage = zeros(Nt,Nblk);
blk.coh = nan(Nt,Nblk);
blk.Surface = nan(1,Nblk);
blk.blend_fringes = zeros(1,Nblk);

do_block_fringe = opts.blend_coreg_en && info.phase_is_unwrapped ...
  && ~isempty(info.dtau_coreg);

for b = 1:Nblk
  cols = blk.starts(b):min(blk.starts(b)+opts.block_size-1,Nx);
  blk.cols{b} = cols;
  wsum = sum(w(:,cols),2);
  blk.dtau(:,b) = sum(w(:,cols).*dtau_w(:,cols),2) ./ wsum;
  blk.coverage(:,b) = sum(info.coh_mask(:,cols),2) / numel(cols);
  blk.coh(:,b) = mean(map.coherence(:,cols),2);
  blk.Surface(b) = mean(map.Surface(cols),'omitnan');

  if do_block_fringe
    diff_blk = info.dtau_coreg(:,cols) - dtau(:,cols);
    sig = info.coh_mask(:,cols) & isfinite(diff_blk);
    if nnz(sig) > 100
      blk.blend_fringes(b) = round(median(diff_blk(sig))*map.fc);
      blk.dtau(:,b) = blk.dtau(:,b) + blk.blend_fringes(b)/map.fc;
    end
  end
end

end
