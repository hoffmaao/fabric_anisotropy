function [dlam, out] = invertHorizontalFabricJoint(obs, par, opts)
%INVERTHORIZONTALFABRICJOINT Regularized joint horizontal-fabric inversion.
%   [dlam, out] = INVERTHORIZONTALFABRICJOINT(obs, par, opts) solves for
%   the piecewise-constant horizontal fabric contrast dlam(z) = lam_x -
%   lam_y over all reflector intervals SIMULTANEOUSLY, minimizing
%       sum_k w_k (dtau_model_k - dtau_obs_k)^2 + alpha ||D dlam||^2
%   where D takes first differences between adjacent intervals. Unlike the
%   exact layer stripping of ptt.invertHorizontalFabric, the smoothness
%   penalty suppresses the noise-amplified oscillation (and eigenvalue
%   bound-pegging) that per-interval exact solves suffer when the
%   traveltime increment across an interval is below the noise level, at
%   the cost of some depth resolution.
%
%   obs fields as in ptt.invertHorizontalFabric (L, z, dtau [ns]), plus
%   optional obs.w: relative misfit weights per node (e.g. coherence);
%   defaults to equal weights.
%   par: firn-ice column model (ptt.defaultParams).
%   opts fields (optional):
%     reg (0.05)  regularization strength; alpha = reg*trace(A'WA)/m so it
%                 is dimensionless relative to the mean data sensitivity
%
%   Outputs match ptt.invertHorizontalFabric: dlam per interval (shallowest
%   first), out.ztop/zbot [m above bed], out.dtau_fit [ns], plus out.alpha,
%   out.rms [ns], out.clipped (intervals clipped to the eigenvalue bound).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'reg') || isempty(opts.reg)
  opts.reg = 0.05;
end

[zs, order] = sort(obs.z(:), 'descend'); % shallowest reflector first
dts = obs.dtau(:);
dts = dts(order);
m = numel(zs);
assert(all(zs > 0 & zs < par.H), 'Reflector heights must satisfy 0 < z < H.');

if isfield(obs,'w') && ~isempty(obs.w)
  w = obs.w(:); w = w(order);
  w(~isfinite(w) | w < 0) = 0;
  if ~any(w > 0), w = ones(m,1); end
else
  w = ones(m,1);
end
W = diag(w);

zhat_nodes = [zs(end:-1:1) / par.H; 1];

% Assumed lam_z at interval midpoints (step approximation), as in the
% layer-stripping variant
zhat_mid = [(1 + zs(1)/par.H)/2; (zs(1:end-1) + zs(2:end)) / (2*par.H)];
lam_z_mid = par.lam_z_bed + (par.lam_z_sfc - par.lam_z_bed) .* zhat_mid;
h_mid = 1 - lam_z_mid;

  function par2 = buildPar(dl)
    dl_nodes = [dl(end:-1:1); dl(1)];
    lz = [lam_z_mid(end:-1:1); lam_z_mid(1)];
    par2 = par;
    par2.fabric = struct('method', 'previous', 'zhat', zhat_nodes, ...
      'lam_x', (([h_mid(end:-1:1); h_mid(1)]) + dl_nodes)/2, 'lam_z', lz);
  end

forward = @(dl) ptt.twttDifference(buildPar(dl), obs.L*ones(m,1), zs);

% Gauss-Newton with Tikhonov first-difference regularization. The forward
% model is nearly linear in dlam, so 2-3 iterations converge.
D = diff(eye(m));
x = zeros(m,1);
delta = 1e-3;
alpha = NaN;
for it = 1:3
  f0 = forward(x);
  if any(~isfinite(f0))
    error('ptt:invertHorizontalFabricJoint:forward', ...
      'Forward model returned non-finite dtau (reflector unreachable?).');
  end
  r = f0 - dts;
  A = zeros(m, m);
  for j = 1:m
    xp = x; xp(j) = xp(j) + delta;
    A(:, j) = (forward(xp) - f0) / delta;
  end
  alpha = opts.reg * trace(A'*W*A) / m;
  dx = -(A'*W*A + alpha*(D'*D)) \ (A'*W*r + alpha*(D'*D)*x);
  x = x + dx;
  % Keep lam_x and lam_y in [0, h]: |dlam| < h
  x = min(max(x, -h_mid + 1e-6), h_mid - 1e-6);
  if norm(dx) < 1e-4
    break;
  end
end

dlam = x;
out.ztop = [par.H; zs(1:end-1)];
out.zbot = zs;
out.dtau_fit = forward(x);
out.alpha = alpha;
out.rms = sqrt(mean(w .* (out.dtau_fit - dts).^2) / mean(w));
out.clipped = abs(abs(x) - h_mid) < 1e-5;
out.par = buildPar(x);

end
