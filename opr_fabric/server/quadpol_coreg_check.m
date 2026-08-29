%SECREG_CHECK Coherence of the ALREADY-coregistered pair in the product.
%
% polarimetric_task.m stores `ref` (HH) and `sec_reg` (VV resampled onto HH
% by the toolbox coregistration), so the co-polarized pair needs no
% recomputation at all - the expensive step has already been paid for and
% shipped. This measures |C_HHVV| from those two directly, per trace and
% with no cross-trace averaging, which is the number the raw-standardphase
% estimate of 0.291 should be compared against.
% repo root derived from this script's own location - see fabric_paths
addpath(fabric_paths());
root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
fn = fullfile(root, 'CSARP_polarimetric', '20250108_02', ...
  'Data_20250108_02_009.mat');
w = {whos('-file', fn).name};
sel = intersect({'ref','sec','sec_reg','Time','Surface', ...
  'row_offset','col_offset'}, w);
P = load(fn, sel{:});

Tv = P.Time(:);
C0 = 299792458; C_ICE = C0/sqrt(3.171);
st = median(P.Surface(:), 'omitnan');
z = (Tv - st)*C_ICE/2;
band = z > 200 & z < 1200;
NRW = 101; kr = ones(NRW,1)/NRW;

A = P.ref; if isstruct(A), A = A.Data; end
fprintf('fields present: %s\n', strjoin(w, ', '));
fprintf('ref size %s\n', mat2str(size(A)));
if isfield(P, 'sec')
  Bs = P.sec; if isstruct(Bs), Bs = Bs.Data; end
  c = H_coh(A, Bs, kr, band);
  fprintf('|C_HHVV| with RAW sec       : median %.3f  p90 %.3f\n', c(1), c(2));
end
if isfield(P, 'sec_reg')
  Br = P.sec_reg;
  c = H_coh(A, Br, kr, band);
  fprintf('|C_HHVV| with COREG sec_reg : median %.3f  p90 %.3f\n', c(1), c(2));
end
if isfield(P, 'row_offset')
  ro = P.row_offset(:);
  fprintf('row_offset median %.3f samples (p5 %.2f p95 %.2f)\n', ...
    median(ro(isfinite(ro))), prctile(ro(isfinite(ro)),5), ...
    prctile(ro(isfinite(ro)),95));
end
if isfield(P, 'col_offset')
  co = P.col_offset(:);
  fprintf('col_offset median %.3f traces\n', median(co(isfinite(co))));
end
fprintf('\nfor comparison, raw standardphase HH/VV uncoregistered: 0.291\n');

function c = H_coh(a, b, kr, band)
num = conv2(a .* conj(b), kr, 'same');
d1 = conv2(abs(a).^2, kr, 'same');
d2 = conv2(abs(b).^2, kr, 'same');
cc = abs(num) ./ sqrt(max(d1.*d2, realmin));
c = [median(cc(band,:), 'all', 'omitnan'), prctile(cc(band,:), 90, 'all')];
end
