%COH_TEST Why is |C_HHVV| only ~0.2, and is it the trace averaging?
%
% Ershadi et al. gate on |C_HHVV| > 0.4 and this survey never reaches it.
% Two candidates, which this separates:
%   (a) the co-polarized channels genuinely decorrelate with along-track
%       DISTANCE - which is what separate H and V antenna phase centres
%       would do, since the two then look at slightly different ice and the
%       common layer phase no longer cancels in the product;
%   (b) the coherence is fine per trace and it is the averaging over 200
%       traces that destroys it.
% If |C| is high per trace and falls as traces are added, it is (a)/(b) and
% the fix is to average less; if it is low even on one trace, the channels
% were never coherent and no amount of care in the averaging helps.
% repo root derived from this script's own location - see fabric_paths
addpath(fabric_paths());
site_root = '/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2';
day_seg = '20250108_02'; frm = 9;
CHAN = {'hh','vv','hv','vh'};
name = sprintf('Data_%s_%03d.mat', day_seg, frm);
S = struct();
for k = 1:4
  q = load(fullfile(site_root, ['CSARP_standardphase_' upper(CHAN{k})], ...
    day_seg, name), 'Data', 'Time', 'Surface', 'Latitude', 'Longitude');
  S.(CHAN{k}) = q.Data;
  if k == 1, ref = q; end
end
Tv = ref.Time(:);
C0 = 299792458; C_ICE = C0/sqrt(3.171);
surf_t = median(ref.Surface(:), 'omitnan');
if ~isfinite(surf_t)
  pw = abs(S.hh).^2; thr = 0.05*max(pw,[],1);
  st = nan(1,size(S.hh,2));
  for j = 1:size(S.hh,2)
    i0 = find(pw(:,j)>thr(j),1); if ~isempty(i0), st(j)=Tv(i0); end
  end
  surf_t = median(st,'omitnan'); clear pw thr;
end
z = (Tv - surf_t)*C_ICE/2;
keep = z>0 & z<1500;
for k=1:4, S.(CHAN{k}) = S.(CHAN{k})(keep,:); end
z = z(keep);
band = z>200 & z<1200;
[~, Nx] = size(S.hh);

% Range-only coherence, averaged over NR range bins, for a growing number
% of traces. NR is held fixed so only the trace count varies.
NR = 101;
kr = ones(NR,1)/NR;
fprintf('frame %s_%03d, %d traces\n', day_seg, frm, Nx);
fprintf('%8s %10s %10s\n', 'ntrace', '|C_hhvv|', '|C_hvvh|');
for nt = [1 2 5 10 25 50 100 200 Nx]
  nt = min(nt, Nx); %#ok<FXSET>
  j = 1:nt;
  num = conv2(sum(S.hh(:,j).*conj(S.vv(:,j)), 2), kr, 'same');
  d1 = conv2(sum(abs(S.hh(:,j)).^2, 2), kr, 'same');
  d2 = conv2(sum(abs(S.vv(:,j)).^2, 2), kr, 'same');
  c = abs(num)./sqrt(max(d1.*d2, realmin));
  numx = conv2(sum(S.hv(:,j).*conj(S.vh(:,j)), 2), kr, 'same');
  e1 = conv2(sum(abs(S.hv(:,j)).^2, 2), kr, 'same');
  e2 = conv2(sum(abs(S.vh(:,j)).^2, 2), kr, 'same');
  cx = abs(numx)./sqrt(max(e1.*e2, realmin));
  fprintf('%8d %10.3f %10.3f\n', nt, median(c(band),'omitnan'), ...
    median(cx(band),'omitnan'));
end

% And against along-track SEPARATION: coherence of trace j with trace j+lag
% in the same channel tells us how fast the ice itself decorrelates, which
% bounds how many traces may legitimately be averaged.
fprintf('\n%8s %12s\n', 'lag', '|C_hh(j,j+lag)|');
for lag = [1 2 5 10 20 50 100]
  if lag >= Nx, break; end
  a = S.hh(:,1:end-lag); b = S.hh(:,1+lag:end);
  num = conv2(sum(a.*conj(b), 2), kr, 'same');
  d1 = conv2(sum(abs(a).^2, 2), kr, 'same');
  d2 = conv2(sum(abs(b).^2, 2), kr, 'same');
  c = abs(num)./sqrt(max(d1.*d2, realmin));
  fprintf('%8d %12.3f\n', lag, median(c(band),'omitnan'));
end
