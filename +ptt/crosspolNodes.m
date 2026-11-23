function N = crosspolNodes(ratio_db, z, opts)
%CROSSPOLNODES The half-wave node ladder in the cross-polarised power 45 deg off the fabric axis.
%
% N = ptt.crosspolNodes(ratio_db, z, opts)
%
% For a birefringent column the cross-polarised power seen 45 deg off the
% principal axes goes as sin^2(delta/2), delta the accumulated two-way
% birefringent phase (Fujita et al. 2006): it is largest where delta is an
% odd multiple of pi (a half-wave plate converts the polarisation
% completely) and vanishes where delta is an even multiple. Each extremum
% of the cross/co power ratio against depth is therefore a node at a KNOWN
% delay, n * (half a carrier fringe), with no phase tracking and no
% unwrapping. On Ridge A the nodes are 150 m apart and the ratio swings
% 2-6 dB. This returns the CANDIDATE extrema; which are nodes, and their
% count n, is decided against the carrier phase by the caller (nodeFabric):
% speckle leaves wiggles of up to a dB or two on the broad top of a
% maximum, and no prominence threshold separates those from a weak node.
%
% Inputs
%   ratio_db  [Nt x 1] cross/co power ratio in dB at the 45-deg azimuth,
%             pooled over traces
%   z         [Nt x 1] depth [m], increasing
%   opts      (all optional)
%     .smooth_m     (20)   running-mean length before extrema are picked
%     .z_min        (200)  no node is accepted above this: the firn and the
%                          near-surface are the ratio's rise from its floor
%     .min_prom_db  (0.5)  an extremum must stand this far from the
%                          neighbouring opposite extrema on both sides
%     .flat_db      (0.15) the node's depth uncertainty is half the span
%                          over which the smoothed ratio stays within this
%                          of the extremum
%
% Output N
%   .z        [n x 1] node depths [m]
%   .zsig     [n x 1] depth uncertainty [m]
%   .type     [n x 1] +1 a maximum (odd multiple of pi), -1 a minimum
%   .prom_db  [n x 1] prominence
%   .smoothed [Nt x 1] the ratio the extrema were read from
if nargin < 3, opts = struct(); end
smooth_m = H_opt(opts, 'smooth_m', 20);
z_min = H_opt(opts, 'z_min', 200);
min_prom = H_opt(opts, 'min_prom_db', 0.5);
flat_db = H_opt(opts, 'flat_db', 0.15);
z = z(:); r = ratio_db(:);
dz = median(diff(z));
k = max(1, round(smooth_m / dz));
ok = isfinite(r);
rs = r;
rs(ok) = movmean(r(ok), k);
rs(~ok) = NaN;

% every local extremum of the smoothed ratio below z_min, alternating
% maxima and minima by construction
d = sign(diff(rs));
d(d == 0) = NaN;
d = fillmissing(d, 'previous');
turn = find(d(1:end-1) .* d(2:end) < 0) + 1;
turn = turn(z(turn) >= z_min & isfinite(rs(turn)) & turn > k & turn <= numel(rs) - k);   % not within a smoothing length of the ends
typ = d(turn - 1);                  % rising before the turn = a maximum (+1)
% pairwise cancellation: speckle leaves hundreds of small wiggles between
% the nodes. Repeatedly remove the ADJACENT pair of turns with the
% smallest swing between them, until every swing is at least min_prom.
% Deleting a neighbouring max-min pair keeps the sequence alternating,
% so what survives is the ladder with its prominence guaranteed on both
% sides; a lone first or last turn that only has one side is judged on
% that side.
idx = turn(:); ty = typ(:);
while numel(idx) >= 2
  sw = abs(diff(rs(idx)));
  [smin, i] = min(sw);
  if smin >= min_prom, break; end
  idx([i, i + 1]) = []; ty([i, i + 1]) = [];
end
% prominence: the smaller swing to the neighbouring turns, the record's
% ends standing in where there is no neighbour - so a curve that merely
% rises to the bottom of the record has no node there
i0 = find(z >= z_min & isfinite(rs), 1);
i1 = find(isfinite(rs), 1, 'last');
pr = nan(size(idx));
for i = 1:numel(idx)
  if i > 1, s1 = abs(rs(idx(i)) - rs(idx(i - 1))); else, s1 = abs(rs(idx(i)) - rs(i0)); end
  if i < numel(idx), s2 = abs(rs(idx(i)) - rs(idx(i + 1))); else, s2 = abs(rs(idx(i)) - rs(i1)); end
  pr(i) = min(s1, s2);
end
keep = pr >= min_prom;
idx = idx(keep); ty = ty(keep); pr = pr(keep);
% depth uncertainty: half the span within flat_db of the extremum
zs = nan(size(idx));
for i = 1:numel(idx)
  v = rs(idx(i));
  if ty(i) > 0, in = rs >= v - flat_db; else, in = rs <= v + flat_db; end
  a = idx(i); b = idx(i);
  while a > 1 && in(a - 1), a = a - 1; end
  while b < numel(rs) && in(b + 1), b = b + 1; end
  zs(i) = max(dz, (z(b) - z(a)) / 2);
end
N = struct('z', z(idx), 'zsig', zs, 'type', ty, 'prom_db', pr, 'smoothed', rs);
end

function v = H_opt(o, name, default)
if isfield(o, name) && ~isempty(o.(name)), v = o.(name); else, v = default; end
end
