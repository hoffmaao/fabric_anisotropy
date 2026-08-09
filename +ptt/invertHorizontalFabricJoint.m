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
%   defaults to equal weights. Optional obs.zref [m above bed]: the height
%   of the surface-reference depth the dtau observations were zeroed at.
%   When present, the forward model is DIFFERENCED against the same
%   reference (model increments zref -> node, so whatever fabric sits
%   above zref cancels exactly instead of being fitted to unobserved ice),
%   the shallowest interval's top is reported at zref rather than at the
%   surface, and a CONSTANT-OFFSET nuisance is estimated jointly: any
%   error in the reference value (surface sidelobes, a coregistration
%   step) shifts every node by the same amount, and with a
%   surface-anchored forward that constant is absorbed almost entirely by
%   the shallowest interval's dlam - 0.3 ns of reference error
%   manufactures dlam ~ 0.06 over a 74 m interval, which is exactly the
%   spurious near-surface fabric the fringe check caught at 57 m. The
%   offset and dlam(1) are degenerate up to the priors, so out flags the
%   first interval as reference-degenerate; quote fabric from interval 2
%   down.
%   par: firn-ice column model (ptt.defaultParams).
%   opts fields (optional):
%     reg (0.05)  regularization strength; alpha = reg*trace(A'WA)/m so it
%                 is dimensionless relative to the mean data sensitivity
%     ref_sigma_ns (3.0)  prior std [ns] of the reference-offset nuisance
%                 (only used when obs.zref is given; loose on purpose -
%                 see the default's comment)
%
%   Outputs match ptt.invertHorizontalFabric: dlam per interval (shallowest
%   first), out.ztop/zbot [m above bed], out.dtau_fit [ns], plus out.alpha,
%   out.rms [ns], out.clipped (intervals clipped to the eigenvalue bound),
%   and with obs.zref: out.ref_offset [ns] and out.ref_degenerate
%   (m x 1 logical, true for the interval that shares its information
%   with the offset).

if nargin < 3, opts = struct(); end
if ~isfield(opts,'reg') || isempty(opts.reg)
  opts.reg = 0.05;
end
if ~isfield(opts,'ref_sigma_ns') || isempty(opts.ref_sigma_ns)
  % Deliberately near-free: interval 1 sits above every node, so its dlam
  % and the offset shift the data identically and only the priors split
  % them - the constant goes to whichever is CHEAPER. Routing 0.3 ns via
  % dlam(1) costs alpha*(0.3/4.7)^2 ~ 0.06 at the production alpha, so
  % the offset prior must sit well below that or the reference error
  % lands back in the fabric (measured: sigma 1.0 still leaked ~75% of
  % it). Nothing REAL is at risk from a near-free constant: fabric
  % cannot produce one across nodes except through interval 1, which is
  % exactly the dof declared unmeasurable by ref_degenerate.
  opts.ref_sigma_ns = 3.0;
end
use_ref = isfield(obs,'zref') && ~isempty(obs.zref) && isfinite(obs.zref);

[zs, order] = sort(obs.z(:), 'descend'); % shallowest reflector first
dts = obs.dtau(:);
dts = dts(order);
m = numel(zs);
assert(all(zs > 0 & zs < par.H), 'Reflector heights must satisfy 0 < z < H.');
if use_ref
  zref = obs.zref;
  assert(zref > zs(1) && zref < par.H, ...
    'Reference height must sit above the shallowest reflector.');
end

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

if use_ref
  % Difference the forward model against the same reference depth the
  % observations were zeroed at: model increments zref -> node, so the
  % fabric above zref (which the data cannot see) cancels exactly.
  forward = @(dl) H_diff(ptt.twttDifference(buildPar(dl), ...
    obs.L*ones(m+1,1), [zref; zs]));
else
  forward = @(dl) ptt.twttDifference(buildPar(dl), obs.L*ones(m,1), zs);
end

% Gauss-Newton with Tikhonov first-difference regularization. The forward
% model is nearly linear in dlam, so 2-3 iterations converge. With
% obs.zref the parameter vector carries the reference-offset nuisance as
% its last entry, priored to N(0, ref_sigma_ns^2): the offset is what a
% reference error looks like in the data, and giving it its own dof is
% what keeps it out of the shallowest interval's dlam.
D = diff(eye(m));
x = zeros(m,1);
c = 0;
delta = 1e-3;
alpha = NaN;
for it = 1:3
  f0 = forward(x);
  if any(~isfinite(f0))
    error('ptt:invertHorizontalFabricJoint:forward', ...
      'Forward model returned non-finite dtau (reflector unreachable?).');
  end
  r = f0 + c - dts;
  A = zeros(m, m);
  for j = 1:m
    xp = x; xp(j) = xp(j) + delta;
    A(:, j) = (forward(xp) - f0) / delta;
  end
  if any(~isfinite(A(:)))
    error('ptt:invertHorizontalFabricJoint:jacobian', ...
      'Perturbed forward model returned non-finite dtau (reflector unreachable?); Jacobian is not finite.');
  end
  alpha = opts.reg * trace(A'*W*A) / m;
  if use_ref
    Aa = [A, ones(m,1)];
    R = blkdiag(alpha*(D'*D), 1/opts.ref_sigma_ns^2);
    dxc = -(Aa'*W*Aa + R) \ (Aa'*W*r + R*[x; c]);
    dx = dxc(1:m);
    c = c + dxc(end);
  else
    dx = -(A'*W*A + alpha*(D'*D)) \ (A'*W*r + alpha*(D'*D)*x);
  end
  if any(~isfinite(dx)) || ~isfinite(c)
    error('ptt:invertHorizontalFabricJoint:solve', ...
      'Gauss-Newton normal-equation solve returned a non-finite step (ill-conditioned system).');
  end
  x = x + dx;
  % Keep lam_x and lam_y in [0, h]: |dlam| < h
  x = min(max(x, -h_mid + 1e-6), h_mid - 1e-6);
  if norm(dx) < 1e-4
    break;
  end
end

dlam = x;
if use_ref
  out.ztop = [zref; zs(1:end-1)];
else
  out.ztop = [par.H; zs(1:end-1)];
end
out.zbot = zs;
out.dtau_fit = forward(x) + c;
out.alpha = alpha;
out.rms = sqrt(mean(w .* (out.dtau_fit - dts).^2) / mean(w));
out.clipped = abs(abs(x) - h_mid) < 1e-5;
out.par = buildPar(x);
out.ref_offset = c;
out.ref_degenerate = false(m,1);
if use_ref
  % The offset and the shallowest interval share their information: node 1
  % is the only observation either touches alone. Flag it; quote fabric
  % from interval 2 down.
  out.ref_degenerate(1) = true;
end

end

function d = H_diff(f)
d = f(2:end) - f(1);
end
