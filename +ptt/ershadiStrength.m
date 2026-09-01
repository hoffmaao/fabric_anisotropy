function out = ershadiStrength(fr, z, geom, opts)
%ERSHADISTRENGTH Infer dlam(z) GIVEN the fitted orientation and scattering.
%
% out = ptt.ershadiStrength(fr, z, geom, opts)
%
% PART TWO of the two-part inversion. Part one (ptt.ershadiInverse, their
% Sect. 3.5) infers the GEOMETRY - the fabric orientation theta and the
% reflection ratio r - while ACCEPTING a dlam profile it never revisits.
% This stage closes the loop: with theta and r held at their fitted values,
% it infers the STRENGTH, the horizontal eigenvalue difference dlam(z),
% from the same layered Fujita forward model.
%
% WHY THIS IS A SEPARATE STAGE AND NOT A JOINT FIT. The two halves have
% opposite conditioning, which is exactly why splitting them works:
%   - theta is a GLOBAL quantity in the observables. phi at any row carries
%     the phase of the whole column above it, so an interval's own theta
%     has leverage only in proportion to its share of that phase - below
%     ~20% the stage absorbs upstream error instead of fitting its layer,
%     and MORE data makes it worse (measured in test_ershadi_r).
%   - dlam is LOCAL. It enters as the phase GRADIENT, d(delta)/dz, so an
%     interval's own dlam is what bends the phase across its own rows
%     whatever the stack above contributed. Given theta, each interval is a
%     well-posed 1-D problem, solved top-down so the accumulated phase
%     entering it is already known.
% Fitting both at once would let a wrong theta be paid for with a wrong
% dlam, which is what "accept dlam, then fit theta" quietly does in reverse.
%
% WEIGHTED, NOT GATED. The direct chain's dlam is NaN wherever the
% coherence falls below its |C| > 0.4 quality gate (their E5). At Ridge A,
% median |C| = 0.47, that voids 61% of deep intervals - and a profile that
% is absent cannot be inferred from. Here the misfit is on the COMPLEX
% coherence weighted by |C|^2, so a low-coherence row contributes little
% instead of nothing, and the estimate degrades smoothly rather than
% falling off a cliff. fr.C_raw (ungated) is used for this reason;
% fr.phi is the gated product and is NOT used.
%
% Inputs
%   fr    ptt.ershadiFabric output (needs psi, C_raw, Cmag, and the
%         recorded constants fc/eps_perp/deps/win_m)
%   z     depth axis (m)
%   geom  part-one geometry, i.e. a ptt.ershadiInverse output or any struct
%         carrying .edges, .theta_int and .r_db_int
%   opts  .dlam_max (0.5)      search ceiling, as their bound
%         .n_grid (61)         coarse grid nodes per interval
%         .w_min (1e-3)        minimum summed weight to attempt an interval
%         .fc/.eps_perp/.deps/.win_m   default from fr, as ershadiInverse
%
% Output
%   .dlam_int    [Nint x 1] inferred strength per interval, NaN where the
%                interval had too little weight to support one
%   .dlam_z      interpolated onto z (nearest interval, NaN outside)
%   .w_int       summed coherence weight behind each interval - the number
%                to check before quoting a value
%   .J_int       final misfit per interval
%   .delta_top   accumulated two-way phase entering each interval, the
%                quantity whose growth makes deep intervals easier, not
%                harder, to constrain
%
% See also ptt.ershadiInverse, ptt.fujitaModel, ptt.ershadiFabric.

if nargin < 4, opts = struct(); end
dlam_max = H_opt(opts, 'dlam_max', 0.5);
n_grid = H_opt(opts, 'n_grid', 61);
w_min = H_opt(opts, 'w_min', 1e-3);
% Intervals whose dlam is KNOWN and should be held rather than fitted:
% [Nint x 1], NaN where the interval is to be fitted. Two uses - carrying a
% value established elsewhere (an ice core, a previous pass), and testing
% how far an upstream error propagates, which is the property that decides
% whether staging the strength stage after the geometry is safe at all.
dlam_fixed = H_opt(opts, 'dlam_fixed', []);

% constants come from the chain that produced the observables, never
% restated independently - see ershadiInverse's H_const note
fwd = struct( ...
  'fc', H_const(opts, fr, 'fc', 750e6), ...
  'eps_perp', H_const(opts, fr, 'eps_perp', 3.15), ...
  'deps', H_const(opts, fr, 'deps', 0.034), ...
  'win_m', H_const(opts, fr, 'win_m', 30));

z = z(:);
if ~isfield(fr, 'C_raw')
  error('ptt:ershadiStrength:noCraw', ...
    ['fr has no C_raw - it must come from a ptt.ershadiFabric that ' ...
    'exports the ungated coherence. The gated fr.phi is not usable ' ...
    'here: it is empty exactly where this stage is needed.']);
end
edges = geom.edges(:);
Nint = numel(edges) - 1;
th = geom.theta_int(:);
rdb = geom.r_db_int(:);
if numel(th) ~= Nint || numel(rdb) ~= Nint
  error('ptt:ershadiStrength:geom', ...
    'geom has %d intervals but %d theta / %d r values', ...
    Nint, numel(th), numel(rdb));
end

psi = fr.psi(:).';
C_obs = fr.C_raw;
W = fr.Cmag.^2;                          % |C|^2 weighting, not a threshold
W(~isfinite(W) | ~isfinite(C_obs)) = 0;
C_obs(~isfinite(C_obs)) = 0;

kz = discretize(z, edges);
dlam_int = nan(Nint, 1);
w_int = zeros(Nint, 1);
J_int = nan(Nint, 1);
delta_top = zeros(Nint, 1);
grid_dl = linspace(0, dlam_max, n_grid);

gpd1 = ptt.birefringentPhaseRate(fwd.fc, fwd.eps_perp, fwd.deps);

% --- top-down: interval k is fitted with the intervals above it already
% resolved, so the phase entering it is known and only its own gradient is
% free. Intervals that cannot be resolved stay NaN and, for the purpose of
% propagating phase to the intervals below, contribute their grid-best
% value with the weight that produced it recorded in w_int - a caller can
% see exactly which links in the chain were weak.
if ~isempty(dlam_fixed)
  dlam_fixed = dlam_fixed(:);
  if numel(dlam_fixed) ~= Nint
    error('ptt:ershadiStrength:dlamFixed', ...
      'dlam_fixed has %d values for %d intervals', ...
      numel(dlam_fixed), Nint);
  end
end
is_fixed = false(Nint, 1);
dl_run = zeros(Nint, 1);
for k = 1:Nint
  if ~isempty(dlam_fixed) && isfinite(dlam_fixed(k))
    dl_run(k) = dlam_fixed(k);
    dlam_int(k) = dlam_fixed(k);
    is_fixed(k) = true;
    delta_top(k) = 2 * gpd1 * sum(dl_run(1:k-1) .* diff(edges(1:k)));
    continue                             % held, not fitted
  end
  rows = find(kz == k);
  if isempty(rows)
    continue
  end
  wsum = sum(W(rows, :), 'all');
  w_int(k) = wsum;
  delta_top(k) = 2 * gpd1 * sum(dl_run(1:k-1) .* diff(edges(1:k)));
  if wsum <= w_min
    continue                             % nothing to infer from
  end
  Jg = inf(size(grid_dl));
  for gi = 1:numel(grid_dl)
    dv = dl_run; dv(k) = grid_dl(gi);
    Jg(gi) = H_cost(dv, th, rdb, edges, z(rows), psi, fwd, ...
      C_obs(rows, :), W(rows, :));
  end
  [~, ib] = min(Jg);
  lo = max(0, grid_dl(max(ib-1, 1)));
  hi = min(dlam_max, grid_dl(min(ib+1, numel(grid_dl))));
  f = @(v) H_cost(H_place(dl_run, k, v), th, rdb, edges, z(rows), psi, ...
    fwd, C_obs(rows, :), W(rows, :));
  if hi > lo
    vb = fminbnd(f, lo, hi, optimset('Display', 'off', 'TolX', 1e-5));
  else
    vb = grid_dl(ib);
  end
  dl_run(k) = vb;
  dlam_int(k) = vb;
  J_int(k) = f(vb);
end

dlam_z = nan(numel(z), 1);
for k = 1:Nint
  m = kz == k;
  dlam_z(m) = dlam_int(k);
end

out = struct('dlam_int', dlam_int, 'dlam_z', dlam_z, 'w_int', w_int, ...
  'is_fixed', is_fixed, ...
  'J_int', J_int, 'delta_top', delta_top, 'edges', edges, ...
  'z_int', 0.5*(edges(1:end-1) + edges(2:end)), 'gpd1', gpd1);

end

% ---------------------------------------------------------------- helpers
function J = H_cost(dl_all, th, rdb, edges, zr, psi, fwd, Cobs, Wt)
%H_COST Weighted complex-coherence misfit over interval k's rows.
% The model is evaluated on the full stack so the phase entering these rows
% is right; only interval k's dlam is being varied by the caller.
NL = numel(th);
layers = struct('top_m', num2cell(edges(1:NL).'), ...
  'dlam', num2cell(dl_all(:).'), 'theta', num2cell(th(:).'), ...
  'r_db', num2cell(rdb(:).'));
mo = ptt.fujitaModel(layers, zr, psi, fwd);
d = mo.C - Cobs;
ok = isfinite(d) & Wt > 0;
if ~any(ok(:)), J = inf; return; end
J = sum(Wt(ok) .* abs(d(ok)).^2) / sum(Wt(ok));
end

function v = H_place(v, k, x)
v(k) = x;
end

function v = H_const(o, fr, f, d)
%H_CONST Prefer the caller's override, else what the chain recorded, else
% the documented default - and refuse a silent disagreement between the
% first two, since dlam here must be inferred under the same constants the
% observables were built with.
have_fr = isstruct(fr) && isfield(fr, f) && ~isempty(fr.(f));
if isstruct(o) && isfield(o, f) && ~isempty(o.(f))
  v = o.(f);
  if have_fr && ~isequal(v, fr.(f))
    error('ptt:ershadiStrength:constMismatch', ...
      ['opts.%s = %g disagrees with the %g recorded by ptt.ershadiFabric. ' ...
      'The observables were built under the recorded value; overriding it ' ...
      'here rescales every modelled phase with nothing in the misfit to ' ...
      'reveal it.'], f, v, fr.(f));
  end
elseif have_fr
  v = fr.(f);
else
  v = d;
end
end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
