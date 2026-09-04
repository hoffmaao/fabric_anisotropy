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
%   .dlam_fwd    the dlam the forward model actually PROPAGATED at every
%                interval: the measured value where dlam_int has one, and
%                an interpolated stand-in where the interval abstained.
%                Read it only to see what was assumed in the gaps - a
%                value here that dlam_int reports as NaN is not a
%                measurement (see the top-down block below)
%   .w_int       summed coherence weight behind each interval - the number
%                to check before quoting a value; 0 with .n_rows_int 0
%                means the interval caught no depth rows at all, which is
%                a different failure from zero coherence over rows it did
%   .n_rows_int  rows of z falling in each interval
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
grid_dl = linspace(0, dlam_max, n_grid);

gpd1 = ptt.birefringentPhaseRate(fwd.fc, fwd.eps_perp, fwd.deps);

if ~isempty(dlam_fixed)
  dlam_fixed = dlam_fixed(:);
  if numel(dlam_fixed) ~= Nint
    error('ptt:ershadiStrength:dlamFixed', ...
      'dlam_fixed has %d values for %d intervals', ...
      numel(dlam_fixed), Nint);
  end
end

% --- per-interval row sets, PADDED by the model's coherence window so the
% modelled C over an interval's rows is formed with the same real
% neighbourhood the observed C_raw had, then trimmed back. ptt.fujitaModel
% builds eq. (7) by conv2 with a win_m kernel and 'same' (zero) padding, so
% evaluating on an interval's rows alone damps the model at both ends: at
% the paper's 50 m intervals and the default win_m = 30 m that is most of
% the interval, and the fitted dlam absorbs the bias. ershadiInverse pads
% for exactly this reason; so does this stage now.
n_rows_int = zeros(Nint, 1);
rows_of = cell(Nint, 1); rows_pad = cell(Nint, 1); keep_of = cell(Nint, 1);
dz_row = median(abs(diff(z)));
npad = max(2, ceil(fwd.win_m / max(dz_row, eps)));
for k = 1:Nint
  r_ = find(kz == k);
  rows_of{k} = r_;
  n_rows_int(k) = numel(r_);
  if isempty(r_), rows_pad{k} = r_; keep_of{k} = []; continue; end
  p0 = max(1, r_(1) - npad); p1 = min(numel(z), r_(end) + npad);
  rows_pad{k} = (p0:p1).';
  keep_of{k} = r_ - p0 + 1;
end

% --- top-down: interval k is fitted with the intervals above it already
% resolved, so the phase entering it is known and only its own gradient is
% free.
%
% TWO ARRAYS, as in ershadiInverse, because propagating and reporting are
% different claims. dlam_int is what was MEASURED and keeps its NaN
% wherever an interval had no rows or too little weight. dl_run is what the
% forward model needs in order to accumulate phase past that interval and
% must carry a number for every one of them - and that number is NOT zero.
% Zero is the specific physical claim that the ice there is isotropic; a
% single low-weight interval asserting it leaves delta_top and every deeper
% interval's phase origin short by 2*gpd1*dlam*L, so every deeper dlam is
% then fitted against a wrong origin. The abstaining intervals are filled
% instead by interpolating the ones that WERE resolved (held at the ends,
% never extrapolated), and the column is swept TWICE so the intervals below
% an abstention are fitted against that filled origin rather than against a
% fabricated zero. What was assumed is returned separately as .dlam_fwd.
zc_int = 0.5*(edges(1:end-1) + edges(2:end));
dl_fill = nan(Nint, 1);
is_fixed = false(Nint, 1);
w_int = nan(Nint, 1);
for pass = 1:2
  dlam_int = nan(Nint, 1);
  J_int = nan(Nint, 1);
  delta_top = nan(Nint, 1);
  dl_run = nan(Nint, 1);
  for k = 1:Nint
    delta_top(k) = 2 * gpd1 * sum(dl_run(1:k-1) .* diff(edges(1:k)));
    rows = rows_of{k};
    wsum = 0;
    if ~isempty(rows), wsum = sum(W(rows, :), 'all'); end
    w_int(k) = wsum;                     % reported for held intervals too
    if ~isempty(dlam_fixed) && isfinite(dlam_fixed(k))
      dl_run(k) = dlam_fixed(k);
      dlam_int(k) = dlam_fixed(k);
      is_fixed(k) = true;
      continue                           % held, not fitted
    end
    if isempty(rows) || wsum <= w_min
      dl_run(k) = H_carry(dl_fill, dl_run, k);
      continue                           % nothing to infer from
    end
    f = @(v) H_cost(H_layers(dl_run, dl_fill, k, v), th, rdb, edges, ...
      z(rows_pad{k}), psi, fwd, C_obs(rows, :), W(rows, :), keep_of{k});
    Jg = arrayfun(f, grid_dl);
    [~, ib] = min(Jg);
    lo = max(0, grid_dl(max(ib-1, 1)));
    hi = min(dlam_max, grid_dl(min(ib+1, numel(grid_dl))));
    if hi > lo
      vb = fminbnd(f, lo, hi, optimset('Display', 'off', 'TolX', 1e-5));
    else
      vb = grid_dl(ib);
    end
    dl_run(k) = vb;
    dlam_int(k) = vb;
    J_int(k) = f(vb);
  end
  ok_k = isfinite(dlam_int);
  if ~any(ok_k), break; end                    % nothing resolved to fill from
  fill = H_fill(zc_int, dlam_int, ok_k, dlam_max);
  if isequal(fill, dl_fill), break; end        % the sweep is stationary
  dl_fill = fill;
end
dlam_fwd = dl_run;

dlam_z = nan(numel(z), 1);
for k = 1:Nint
  m = kz == k;
  dlam_z(m) = dlam_int(k);
end

out = struct('dlam_int', dlam_int, 'dlam_z', dlam_z, 'dlam_fwd', dlam_fwd, ...
  'w_int', w_int, 'n_rows_int', n_rows_int, ...
  'is_fixed', is_fixed, ...
  'J_int', J_int, 'delta_top', delta_top, 'edges', edges, ...
  'z_int', zc_int, 'gpd1', gpd1);

end

% ---------------------------------------------------------------- helpers
function J = H_cost(dl_all, th, rdb, edges, zr, psi, fwd, Cobs, Wt, keep)
%H_COST Weighted complex-coherence misfit over interval k's rows.
% The model is evaluated on the full stack so the phase entering these rows
% is right; only interval k's dlam is being varied by the caller. `zr` is
% the interval's rows PADDED by the coherence window, and `keep` trims the
% model back to the interval's own rows after evaluation, which is what
% Cobs and Wt were sliced to; [] keeps all.
NL = numel(th);
layers = struct('top_m', num2cell(edges(1:NL).'), ...
  'dlam', num2cell(dl_all(:).'), 'theta', num2cell(th(:).'), ...
  'r_db', num2cell(rdb(:).'));
mo = ptt.fujitaModel(layers, zr, psi, fwd);
Cm = mo.C;
if nargin >= 10 && ~isempty(keep), Cm = Cm(keep, :); end
d = Cm - Cobs;
ok = isfinite(d) & Wt > 0;
if ~any(ok(:)), J = inf; return; end
J = sum(Wt(ok) .* abs(d(ok)).^2) / sum(Wt(ok));
end

function v = H_layers(dl_run, dl_fill, k, x)
%H_LAYERS The full-column dlam the forward model needs while interval k is
% being fitted. Above k: what the sweep has already resolved or carried.
% At k: the trial value. BELOW k: nothing has been fitted yet, but the
% downward half of the padded row set reaches into interval k+1 and the
% coherence window there is what makes the bottom of interval k's own rows
% modelled correctly. A NaN would not merely be unknown, it would silently
% DROP those rows from the misfit (H_cost skips non-finite cells), so the
% pad would buy nothing at the bottom edge. The previous sweep's
% interpolated profile is used where there is one, else the trial value is
% continued downward - and the second sweep, which has the real profile,
% is what makes this assumption stop mattering.
v = dl_run(:);
v(k) = x;
below = k+1:numel(v);
v(below) = dl_fill(below);
v(below(~isfinite(v(below)))) = x;
end

function v = H_carry(dl_fill, dl_run, k)
%H_CARRY The value an ABSTAINING interval contributes to the phase
% propagated below it. The interpolated fill from the previous sweep if
% there is one, else the nearest resolved interval above; only a column in
% which nothing at all has been resolved yet falls back to 0, and every
% interval of such a column abstains anyway.
if k <= numel(dl_fill) && isfinite(dl_fill(k)), v = dl_fill(k); return; end
j = find(isfinite(dl_run(1:k-1)), 1, 'last');
if isempty(j), v = 0; else, v = dl_run(j); end
end

function fill = H_fill(zc, dlam_int, ok_k, dlam_max)
%H_FILL Resolved intervals interpolated onto every interval centre, HELD at
% the ends rather than extrapolated - a linear extrapolation of dlam past
% the last measured interval runs negative, or past the search ceiling,
% with nothing in the data to stop it.
kk = find(ok_k);
if isscalar(kk)
  fill = repmat(dlam_int(kk), numel(zc), 1);
else
  fill = interp1(zc(kk), dlam_int(kk), zc, 'linear');
  fill(zc < zc(kk(1))) = dlam_int(kk(1));
  fill(zc > zc(kk(end))) = dlam_int(kk(end));
end
fill = min(max(fill, 0), dlam_max);
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
