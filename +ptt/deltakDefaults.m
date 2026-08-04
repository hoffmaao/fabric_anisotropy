function dk = deltakDefaults(opts)
%DELTAKDEFAULTS Defaulted delta-k ladder options (opts.deltak).
%   dk = DELTAKDEFAULTS(opts) returns opts.deltak with every field the
%   split-spectrum ladder uses filled in (see ptt.deltakTraveltime for what
%   each one means). It is the single definition of the analysis-cell
%   geometry: ptt.deltakTraveltime smooths and interpolates on that grid,
%   and ptt.imgCombSeam has to widen the seam band by the same reach, so
%   neither may carry its own copy of the numbers.

if nargin < 1 || ~isstruct(opts) || ~isfield(opts,'deltak') || isempty(opts.deltak)
  dk = struct();
else
  dk = opts.deltak;
end
if ~isfield(dk,'n_sub') || isempty(dk.n_sub), dk.n_sub = 12; end
if ~isfield(dk,'cell_twtt') || isempty(dk.cell_twtt), dk.cell_twtt = 100e-9; end
if ~isfield(dk,'cell_ntr') || isempty(dk.cell_ntr), dk.cell_ntr = 31; end
if ~isfield(dk,'smooth') || isempty(dk.smooth), dk.smooth = 5; end
if ~isfield(dk,'band_frac') || isempty(dk.band_frac), dk.band_frac = [0.02 0.98]; end
if ~isfield(dk,'orientation') || isempty(dk.orientation), dk.orientation = 0; end
if ~isfield(dk,'ref_band') || isempty(dk.ref_band), dk.ref_band = [150e-9 600e-9]; end

end
