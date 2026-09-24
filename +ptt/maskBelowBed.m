function T = maskBelowBed(T, z, z_bed, margin)
%MASKBELOWBED Blank every channel below each trace's own bed.
%
% T = ptt.maskBelowBed(T, z, z_bed, margin)
%
% T       struct of complex [Nt x Nx] channels (hh, vv, hv, vh, or any
%         subset - every field of that shape is masked)
% z       [Nt x 1] depth [m]
% z_bed   [Nx] per-trace bed depth [m]; NaN where the bed is unpicked, and
%         such a trace is left whole (not masking is recoverable, masking
%         ice away is not)
% margin  metres above the bed that are also blanked (default 20), so the
%         bed return's own leading edge does not count as ice
%
% Samples at or below z_bed - margin become NaN. Everything downstream
% averages traces with 'omitnan' (ptt.quadpolMoments), so a moment at depth
% z is then formed only from traces whose ice reaches z, however much the
% bed varies along the line - a single depth cut for a segment cannot do
% that where the ice runs 567-996 m inside one Taylor Dome frame.
%
% WHY. Below the bed the record is not noise: the bed return's tail,
% off-nadir energy and, on the shelves, basal multiples are coherent and
% polarised (measured 23 Sep 2026: HV-VH coherence 0.97 and HH-VV |C|
% 0.5-0.6 up to 150 m below the pick at Taylor Dome; the Eastwind shelf
% base 17 dB above the ice). The fabric model reads that as fabric with an
% axis of its own. See SUB-BED WINDOWS in ptt.quadpolFabricLS.
%
% See also ptt.quadpolFrameTheta, ptt.quadpolMoments.

if nargin < 4 || isempty(margin), margin = 20; end
z = z(:);
fn = fieldnames(T);
Nt = numel(z);
Nx = [];
for k = 1:numel(fn)
  if size(T.(fn{k}), 1) == Nt, Nx = size(T.(fn{k}), 2); break; end
end
if isempty(Nx), return; end
z_bed = z_bed(:).';
if numel(z_bed) ~= Nx
  error('ptt:maskBelowBed:size', 'z_bed has %d entries for %d traces', numel(z_bed), Nx);
end
cut = z_bed - margin;                   % NaN beds give a NaN cut: no mask
mask = z >= cut;                        % [Nt x Nx]; NaN compares false
if ~any(mask(:)), return; end
for k = 1:numel(fn)
  v = T.(fn{k});
  if ~isnumeric(v) || size(v, 1) ~= Nt || size(v, 2) ~= Nx, continue; end
  v(mask) = NaN;
  T.(fn{k}) = v;
end
end
