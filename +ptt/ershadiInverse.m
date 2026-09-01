function out = ershadiInverse(fr, z, opts)
%ERSHADIINVERSE Ershadi et al. (2022) Sect. 3.5 inversion for theta and r.
%
% out = ptt.ershadiInverse(fr, z, opts)
%
% The constrained non-linear least-squares step of Ershadi et al. (2022,
% The Cryosphere 16, 1719-1739): given the observables ptt.ershadiFabric
% extracts - the power anomalies dP_HH and dP_HV (their eq. 12), the
% coherence phase phi_HHVV (eqs. 7-8) - fit the layered Fujita forward
% model (ptt.fujitaModel, their eq. 5) for a piecewise-constant profile of
% the fabric orientation theta_i and the reflection ratio r_i, while the
% horizontal anisotropy profile is ACCEPTED from the phase-gradient
% estimate and not re-fit ("we optimize theta and r for all depth
% intervals, while at this stage we accept the estimated dlam0", 3.5.4).
%
% FAITHFUL CHOICES, with their sources:
%   - Piecewise-constant intervals (their EDML parameterization, 3.5.2;
%     50 m default). The Legendre alternative is not implemented.
%   - Initial guess (3.5.3): theta0 from ptt.ershadiFabric's theta output,
%     which is exactly their recipe (dP_HV minima, disambiguated by the
%     phase polarity and the sign of Psi); r0 = 0 dB.
%   - Cost (3.5.4, eqs. 15-18): J = l1*J_phi + l2*J_dPHH + l3*J_dPHV with
%     each field STANDARDIZED - here by the observation's own mean and
%     standard deviation over finite entries, the same affine map applied
%     to the model, so the three terms are commensurate and a misfit
%     cannot hide in an overall scale.
%   - Weights (Table 3) are 0/1 selectors set PER PARAMETER: theta is fit
%     with w_theta (Dome C row [1 0 0]: coherence phase only; EDML row
%     [0 1 0]: dP_HH only, for strong anisotropy where the phase misfit
%     "was not applicable"), r with w_r (both sites [0 1 0]).
%   - The two stages are CYCLED (theta, then r, repeated) until the total
%     cost stops improving, capped at opts.n_cycles. The paper does not say
%     whether they cycle: Table 3's per-parameter weight rows force the
%     staging, and cycling is the conservative completion of it. theta and
%     r are coupled in both observables - r sets the co-pol node azimuths
%     of dP_HH (that is the whole basis of eq. 13) and it scales the
%     C_HHVV phase excursion - so a theta fitted while r is pinned at the
%     r0 = 0 dB guess is not the theta that goes with the fitted r.
%   - Bounds (3.5.4): 0 < theta_i < pi, -30 dB < r_i < 30 dB, enforced by
%     fmincon bound constraints (they used log-barriers inside the cost
%     with fmincon; interior-point bounds are the same mechanism, owned by
%     the solver instead of hand-rolled).
%   - r = Gamma_y/Gamma_x with r_dB = 20*log10(r) - the amplitude
%     convention their eq. (13) fixes (see ptt.fujitaModel).
%
% Also returned: the eq. (13) ANALYTIC estimate r = 1/tan^2(AD/2). It is
% OPTIMIZER-INDEPENDENT but DATA-CHAIN-CONDITIONED - not "optimization
% free", because it needs an axis to resolve eq. (13)'s arc ambiguity and
% an accumulated phase to find the anti-phase depths. Both are taken from
% ptt.ershadiFabric: the initial-guess theta (th0_int, their 3.5.3 recipe
% of dP_HV minima disambiguated by phase polarity) and the accepted dlam.
% Neither is ever taken from the FITTED theta, and that is the whole point:
%   CAN catch - a theta stage that flips 90 deg or drifts from its initial
%     guess. Conditioning on the fitted theta instead would move the
%     admitted rows and flip the arc in step with the fit, so eq. (13)
%     would return ~1/r exactly where the r stage also returned ~1/r; the
%     cross-check would confirm the error rather than expose it.
%   CANNOT catch - an error already in the initial-guess chain itself. If
%     ershadiFabric's theta or dlam is wrong, r13 is wrong the same way and
%     agrees with a fit that inherited the same wrong start.
% Evaluated only where two co-pol nodes resolve in dP_HH, the depth is near
% ANTI-PHASE, and the stack above is near CO-AXIAL - see the r13 block for
% all three gates and the measurements behind them.
%
% Inputs
%   fr    output struct of ptt.ershadiFabric (needs psi, dP_hh, dP_hv,
%         phi, Cmag, theta, dlam; fc/eps_perp/deps/win_m are read from it
%         when present - see below)
%   z     depth axis (m), same rows as the fields in fr
%   opts  .interval_m (50)   piecewise interval length
%         .z_fit ([200 z(end)]) band the misfits are evaluated over
%         .w_theta ([1 0 0])  Table-3 row for the theta stage
%         .w_r     ([0 1 0])  Table-3 row for the r stage
%         .n_cycles (3)       cap on theta->r staging cycles
%         .cycle_tol (1e-3)   relative total-cost improvement to continue
%         .fit_decim (4)      depth-row decimation inside the cost
%         .fit_psi_decim (4)  azimuth-column decimation inside the cost
%         .fc, .eps_perp, .deps, .win_m   passed to ptt.fujitaModel, and
%              DEFAULTED from the values ptt.ershadiFabric recorded in fr.
%              These constants are already burned into the accepted dlam,
%              so restating one differently in opts is an error, not an
%              override: it would rescale every modelled phase with nothing
%              in the cost to reveal it.
%         .r13_cos_max (-0.85) anti-phase gate for the eq.-(13) estimate
%         .r13_coax_deg (15)  co-axial gate for the eq.-(13) estimate: how
%              far, modulo 90 deg, two interval axes above the row may lie
%              before the scalar phase accumulation stops being meaningful
%         .r13_neg_frac (0.2) an interval is NEGLIGIBLE, and so exempt from
%              both the co-axial test and the phase sum, when its worst-case
%              contribution is under this fraction of the anti-phase
%              half-window; if the exempted total reaches that same bound
%              the exemption is withdrawn and they are checked instead
%         .r13_axis_res_min (0.5) an interval's axis counts as DETERMINED
%              when its doubled-angle resultant length reaches this; below
%              it the axis is unknowable and the column below abstains
%         .max_iter (150)     fmincon iteration cap per stage
%
% Output
%   .z_int, .theta_int, .r_db_int   per-interval fitted values
%   .theta, .r_db                   the same, expanded to the z rows
%   .dlam_int                       the accepted dlam per interval
%   .r13_db                         eq.-(13) analytic r per depth row, NaN
%                                   wherever any gate abstained
%   .r13_reason                     per-row code naming the gate that fired,
%                                   indexing .r13_reason_key: 0 resolved,
%                                   1 coaxial veto, 2 undefined axis,
%                                   3 demoted small intervals unreconciled,
%                                   4 off anti-phase, 5 nodes unresolvable.
%                                   NaN outside the band. 1-3 are properties
%                                   of the column, 4 of the depth, 5 of the
%                                   data - they are not interchangeable
%   .r13_reason_key                 cellstr naming codes 0..5
%   .r13_gate                       the per-interval code 0..3 behind them
%   .theta0_int                     the INITIAL-GUESS theta per interval -
%                                   returned because r13 is conditioned on
%                                   it, so measuring theta_int against it is
%                                   how a caller confirms that an r13-vs-fit
%                                   disagreement really is a drifted or
%                                   flipped theta stage
%   .theta0_res                     doubled-angle resultant length behind
%                                   each theta0_int (0 = no axis at all)
%   .J0, .J                         cost immediately before/after each
%                                   stage, [n_cycles_used x 2], columns
%                                   [theta stage, r stage] - so each pair
%                                   brackets that stage alone
%   .exitflag                       fmincon exit flags, same shape
%   .n_cycles_used                  staging cycles actually run
%   .cycle_best                     index of the cycle whose end point is
%                                   returned; < n_cycles_used means a later
%                                   cycle was worse and was discarded
%   .cycle_worsened                 true if any cycle raised the combined
%                                   staged cost
%
% See also ptt.fujitaModel, ptt.ershadiFabric.

if nargin < 3, opts = struct(); end
int_m = H_opt(opts, 'interval_m', 50);
zfit = H_opt(opts, 'z_fit', [200, z(end)]);
w_th = H_opt(opts, 'w_theta', [1 0 0]);
w_r = H_opt(opts, 'w_r', [0 1 0]);
dec_z = H_opt(opts, 'fit_decim', 4);
dec_p = H_opt(opts, 'fit_psi_decim', 4);
max_iter = H_opt(opts, 'max_iter', 150);
n_cycles = H_opt(opts, 'n_cycles', 3);
cycle_tol = H_opt(opts, 'cycle_tol', 1e-3);
r13_cos_max = H_opt(opts, 'r13_cos_max', -0.85);
% Out of range this fails SILENTLY rather than loudly: acos(-1.5) is
% complex, so the half-window and the negligible bound derived from it are
% complex too, and MATLAB's < then compares real parts only - nothing is
% ever negligible, the closure check never fires, and the anti-phase test
% never abstains, so the gate quietly disables itself with no error.
if ~isscalar(r13_cos_max) || ~isreal(r13_cos_max) || ...
    ~isfinite(r13_cos_max) || r13_cos_max <= -1 || r13_cos_max >= 1
  error('ptt:ershadiInverse:r13CosMax', ...
    'r13_cos_max must be a real scalar in (-1, 1), got %s', ...
    mat2str(r13_cos_max));
end
r13_coax_deg = H_opt(opts, 'r13_coax_deg', 15);
r13_neg_frac = H_opt(opts, 'r13_neg_frac', 0.2);
r13_axis_res_min = H_opt(opts, 'r13_axis_res_min', 0.5);

% The forward constants are NOT free here: fr.dlam was produced by
% ptt.ershadiFabric under a particular fc/eps_perp/deps and is accepted
% rather than re-fit, so the model must propagate phase under the same
% ones. H_const takes them from fr and rejects a disagreeing opts value.
fwd = struct('fc', H_const(opts, fr, 'fc', 750e6), ...
  'eps_perp', H_const(opts, fr, 'eps_perp', 3.15), ...
  'deps', H_const(opts, fr, 'deps', 0.034), ...
  'win_m', H_const(opts, fr, 'win_m', 30));

gpd1 = ptt.birefringentPhaseRate(fwd.fc, fwd.eps_perp, fwd.deps);

z = z(:);
band = z >= zfit(1) & z <= zfit(2);
if ~any(band)
  error('ptt:ershadiInverse:band', 'z_fit selects no rows');
end

% --- intervals over the fit band; the layer stack starts at the SURFACE
% (isotropic-ish shallow rows still propagate phase), so intervals above
% the band exist too, carrying the accepted dlam and the initial theta.
z0 = 0; z1 = min(zfit(2), z(end));
edges = (z0:int_m:z1).';
if edges(end) < z1, edges(end+1) = z1; end
Nint = numel(edges) - 1;
% Exact centres, not edges + int_m/2: the append above can leave a SHORT
% final interval whose centre is not half an interval below its top.
zc_int = 0.5 * (edges(1:end-1) + edges(2:end));
% Single owner of "which interval does this row belong to". discretize
% closes the LAST bin on the right, so the row at exactly edges(end) lands
% in interval Nint instead of falling through every half-open z < edges(k+1)
% test and silently becoming NaN in dlam_int's median, in theta/r_db, and in
% the r13 reason accounting.
kz = discretize(z, edges);
fitset = find(zc_int >= zfit(1));                    % intervals we optimize
if isempty(fitset)
  error('ptt:ershadiInverse:fitset', ...
    ['no interval centre reaches z_fit(1) = %g m (interval_m = %g m over ' ...
     '0-%g m); widen z_fit or shorten interval_m'], zfit(1), int_m, z1);
end

% --- accepted dlam and initial theta per interval (3.5.3/3.5.4)
dlam_int = zeros(Nint, 1);
th0_int = zeros(Nint, 1);
th0_res = zeros(Nint, 1);
for k = 1:Nint
  m = kz == k;
  d = fr.dlam(m); d = d(isfinite(d));
  if ~isempty(d), dlam_int(k) = median(d); end
  t = fr.theta(m); t = t(isfinite(t));
  if ~isempty(t)
    % the doubled-angle resultant LENGTH is the direct measure of whether
    % this interval has an axis at all: ptt.ershadiFabric resolves theta's
    % 90 deg polarity branch row by row from sign(psi_grad), so an interval
    % whose rows split between the branches cancels to a near-zero
    % resultant at an essentially arbitrary angle
    R2 = mean(exp(2i * t));
    th0_int(k) = mod(angle(R2) / 2, pi);
    th0_res(k) = abs(R2);
  elseif k > 1
    th0_int(k) = th0_int(k-1);
  end
end
r0_int = zeros(Nint, 1);                 % "The initial guess for r0dB is zero"

% --- observation fields, decimated and standardized once
iz = find(band); iz = iz(1:dec_z:end);
ip = 1:dec_p:numel(fr.psi);
psi_fit = fr.psi(ip);
obs.phi = fr.phi(iz, ip);
obs.hh = fr.dP_hh(iz, ip);
obs.hv = fr.dP_hv(iz, ip);
zsub = z(iz);
std_of = @(A) deal_std(A);
[mu1, s1] = std_of(obs.phi); [mu2, s2] = std_of(obs.hh); [mu3, s3] = std_of(obs.hv);
nrm = {@(A) (A - mu1)/s1, @(A) (A - mu2)/s2, @(A) (A - mu3)/s3};
obs_n = {nrm{1}(obs.phi), nrm{2}(obs.hh), nrm{3}(obs.hv)};

cost = @(th_all, r_all, w) H_cost(th_all, r_all, dlam_int, edges, ...
  zsub, psi_fit, fwd, obs_n, nrm, w);

oopt = optimoptions('fmincon', 'Display', 'off', 'Algorithm', ...
  'interior-point', 'MaxIterations', max_iter, ...
  'MaxFunctionEvaluations', 200 * numel(fitset));

th = th0_int; rdb = r0_int;
nf = numel(fitset);
lb_th = zeros(nf, 1);      ub_th = pi * ones(nf, 1);      % 0 < theta < pi
lb_r = -30 * ones(nf, 1);  ub_r = 30 * ones(nf, 1);       % +-30 dB

% --- staged fit, cycled. Each J0/J pair brackets its own stage: J0 is
% evaluated at the parameters the stage STARTS from, so J0 - J is that
% stage's improvement and nothing else.
%
% BEST-SO-FAR, not last. A cycle is not a descent step on the combined
% cost: stage 1 minimizes the w_theta cost at the PREVIOUS r, which can
% raise the w_r cost, and stage 2 then only guarantees J(cyc,2) <=
% J0(cyc,2) = the w_r cost at the new theta and the old r - a quantity that
% may already exceed the previous cycle's J(cyc,2). The same holds for the
% w_theta cost across the r update. So a cycle can end strictly worse than
% the one before it, and the break below fires on "did not improve by
% cycle_tol", which is true for a worsening cycle too. Keeping the best
% (th, rdb, Jtot) triple means convergence and divergence do not have to be
% told apart by luck: the worse iterate is discarded either way.
J0 = nan(n_cycles, 2); J = nan(n_cycles, 2); ex = zeros(n_cycles, 2);
Jtot_prev = Inf;
Jtot_best = Inf; th_best = th; rdb_best = rdb; cyc_best = 0;
worsened = false;
n_used = 0;
for cyc = 1:n_cycles
  n_used = cyc;

  % stage 1: theta intervals against w_theta (Table 3, theta row), r held
  if isnan(J0(cyc, 1)), J0(cyc, 1) = cost(th, rdb, w_th); end
  f1 = @(p) cost(H_place(th, fitset, p), rdb, w_th);
  [p1, J(cyc,1), ex(cyc,1)] = fmincon(f1, th(fitset), [], [], [], [], ...
    lb_th, ub_th, [], oopt);
  th = H_place(th, fitset, p1);

  % stage 2: r intervals against w_r (Table 3, r row), theta held
  J0(cyc, 2) = cost(th, rdb, w_r);
  f2 = @(p) cost(th, H_place(rdb, fitset, p), w_r);
  [p2, J(cyc,2), ex(cyc,2)] = fmincon(f2, rdb(fitset), [], [], [], [], ...
    lb_r, ub_r, [], oopt);
  rdb = H_place(rdb, fitset, p2);

  Jth_end = cost(th, rdb, w_th);
  Jtot = Jth_end + J(cyc, 2);
  if Jtot < Jtot_best
    Jtot_best = Jtot; th_best = th; rdb_best = rdb; cyc_best = cyc;
  end
  worsened = worsened || Jtot > Jtot_prev;

  % converged when a whole cycle no longer buys a relative cycle_tol of the
  % summed staged cost at the cycle's own end point.
  % isfinite guard: the first cycle has no predecessor, and Inf - J <= Inf
  % would otherwise be true and break before any cycling happened.
  if isfinite(Jtot_prev) && ...
      Jtot_prev - Jtot <= cycle_tol * max(abs(Jtot_prev), realmin)
    break;
  end
  Jtot_prev = Jtot;
  if cyc < n_cycles, J0(cyc+1, 1) = Jth_end; end
end
th = th_best; rdb = rdb_best;
J0 = J0(1:n_used, :); J = J(1:n_used, :); ex = ex(1:n_used, :);

% --- eq. (13) analytic r from the co-pol node angular distance, where two
% nodes are resolvable in a depth row of dP_HH. The node pair defines two
% arcs (AD and pi - AD) and eq. (13) applied to the wrong one returns 1/r -
% the axis says which arc is which, since the nodes straddle it.
%
% CONDITIONED ON th0_int, NEVER ON THE FITTED theta. Both the arc below and
% the phase accumulation use the initial-guess axis, which ptt.ershadiFabric
% derived from the data alone. Using the fitted theta would make this agree
% with the fit by construction in exactly the case worth catching: theta is
% bounded [0, pi] and the EDML w_theta = [0 1 0] configuration fits it
% against dP_HH, which is near-symmetric under theta -> theta + 90, so the
% stage can land an interval 90 deg off. At a flipped axis every
% cos(2*dtheta) below changes sign, relocating the admitted rows, AND the
% arc test flips so eq. (13) returns 1/r - while the r stage, fitting at
% that same flipped axis, returns ~1/r too. The two errors would cancel and
% verdict the flip as agreement. Off th0_int they do not cancel, so the
% flip shows up as an r13-vs-fit disagreement.
%
% RESTRICTED TO ANTI-PHASE. r = 1/tan^2(AD/2) holds only where the two-way
% birefringent phase delta is pi: there |s_hh|^2 = (cos^2 - r sin^2)^2 has
% true nulls at tan^2 = 1/r. Off anti-phase the minimum of
%   |s_hh|^2 = cos^4 + r^2 sin^4 + 2 r cos^2 sin^2 cos(delta)
% moves, and eq. (13) returns (r^2 - r cos delta)/(1 - r cos delta) instead
% of r. For a +10 dB truth that reads +10.0 dB at cos delta = -1, +13.0 dB
% at -0.5 and +20.0 dB at 0 - and at quadrature the minimum is still ~4.8 dB
% below the row median, so a depth-of-null gate does NOT stand in for a
% phase gate and the estimate is biased high wherever it resolves off
% anti-phase. delta(z) is therefore accumulated from the ACCEPTED dlam
% profile and rows are admitted only where cos(delta) <= r13_cos_max
% (-0.85 by default: |delta - pi| < 32 deg, capping the residual bias near
% 0.7 dB at +10 dB). This is also the faithful reading - Ershadi measure AD
% at the co-pol nodes, which is where the anti-phase depths are.
%
% Accumulating delta needs a COMMON FRAME. ptt.ershadiFabric reports dlam
% as a magnitude (eq. 10 is evaluated at whichever principal axis carries
% lam_max), so an interval whose axis lies 90 deg from the reference one
% has its fast and slow axes swapped and contributes with the OPPOSITE
% sign. Summing magnitudes would place the anti-phase depths wrongly in
% exactly the geometry this retrieval is aimed at - the EDML-shaped profile
% whose deep zone sits at theta + 90 - and eq. (13) would then be read off
% node. Each interval is therefore projected onto the reference interval's
% axis with cos(2*dtheta).
%
% AND THAT PROJECTION IS ONLY VALID NEAR 0 OR 90 DEG, so rows whose stack
% above is not near co-axial ABSTAIN. cos(2*dtheta) is EXACT at 0 and 90
% deg - R(90) diag(e^{ix}, e^{-ix}) R(90)' = diag(e^{-ix}, e^{ix}), a clean
% sign flip - but it is degenerate in between, not merely approximate. At
% dtheta = 45 deg it contributes exactly zero, because
% R(45) diag(e^{ix}, e^{-ix}) R(45)' = [cos x, i sin x; i sin x, cos x]:
% the accumulated DIAGONAL phase difference really is zero, but the layer's
% whole effect has moved into the off-diagonal coupling that this scalar
% accumulation discards, and "anti-phase" stops predicting where the co-pol
% nodes sit at all. Modelling that coupling is not the job of a
% cross-check, so the gate abstains instead: every pair of contributing
% intervals above the row must lie within r13_coax_deg of each other modulo
% 90 deg. Modulo 90 because 90 deg IS co-axial - same eigen-axes, swapped
% labels - which is what keeps the EDML-shaped deep zone admissible.
% CONSEQUENCE, stated plainly: at sites where the axis rotates gradually
% with depth - Dome C, EDML, Thwaites - r13 abstains over most of the
% column. That is correct behaviour, not a defect; the fitted profile is
% the product there and r13 simply has nothing valid to say.
% The gate removes the regime where the scalar accumulation is meaningless.
% It does NOT bound the residual accumulated-phase error, which grows with
% depth and with how much phase the off-axis intervals carry - one more
% reason r13 is a cross-check and not a measurement.
%
% ONE DEFINITION OF NEGLIGIBLE, used by both gates. An interval is
% negligible when its OWN worst-case phase contribution dphi_int is small
% against the quantity the gate protects - the anti-phase half-window
% w_ap = pi - acos(r13_cos_max) - and never as a percentage of the column
% total. A relative rule is what makes the two gates disagree with each
% other: at 1% of this test's own 51.2 rad budget the exemption is 0.51 rad
% = 29.3 deg, against a half-window of 0.555 rad = 31.8 deg, so an interval
% could shift the gate by 92% of its half-width while never having to be
% co-axial with anything.
% A negligible interval is then treated consistently on BOTH sides: it is
% dropped from cum_dl (so it contributes no unchecked phase) and it is
% exempt from the pairwise co-axiality test (an interval that cannot move
% the phase cannot decorrelate the projection either, so it must not veto).
% CLOSURE, which is the whole point: dropping them is only sound while
% their contributions cannot ADD UP, so the SUM of all exempted worst cases
% is held under the same bound. Exceed it and the exemption is WITHDRAWN -
% every small interval is demoted back into the checked set, so many small
% unchecked shifts can never accumulate into one large one.
% The reference interval k is never exempt whatever phase it carries - its
% axis is the frame cum_dl projects onto and the arc H_r13 reads.
%
% NEGLIGIBLE PHASE AND UNDEFINED AXIS ARE DIFFERENT, and only the first is
% an exemption. An interval with meaningful phase but no determined axis
% (th0_res below r13_axis_res_min - its rows split across ershadiFabric's
% 90 deg polarity branches, the weak-fabric case) genuinely has an
% unknowable projection, so it MUST abstain the column below it. That is
% correct abstention, not a defect. The two report separately below.
Lint = diff(edges);
dphi_int = 2 * gpd1 * dlam_int .* Lint;   % two-way phase each interval adds
coax_min = cos(2 * deg2rad(r13_coax_deg));
w_ap = pi - acos(r13_cos_max);            % anti-phase half-window, rad
neg_bound = r13_neg_frac * w_ap;

cum_dl = zeros(Nint, 1);
gate = zeros(Nint, 1);
for k = 1:Nint
  small = dphi_int(1:k) < neg_bound;
  small(k) = false;
  % Closure failure means STOP EXEMPTING, not abstain. Discarding the row
  % would throw away coverage the co-axial test might well have granted -
  % four 50 m intervals at dlam = 0.006 each sit under the bound yet total
  % 0.36 rad, and they may be perfectly co-axial - and it would leave the
  % row's delta2 short by exactly the amount that tripped the gate. Demoting
  % them into the active set restores the full phase AND subjects them to
  % the same axis and co-axiality tests, so nothing unchecked accumulates
  % either way; only a demoted set that then fails those tests abstains.
  demoted = sum(dphi_int(small)) >= neg_bound;
  if demoted, small(:) = false; end
  act = find(~small);
  jj = act(act < k);
  cum_dl(k) = sum(dlam_int(jj) .* Lint(jj) .* ...
    cos(2 * (th0_int(jj) - th0_int(k))));
  ta = th0_int(act);
  ok_axis = all(th0_res(act) >= r13_axis_res_min);
  ok_coax = all(abs(cos(2 * (ta - ta.'))) >= coax_min, 'all');
  if ok_axis && ok_coax
    gate(k) = 0;
  elseif demoted
    gate(k) = 3;
  elseif ~ok_axis
    gate(k) = 2;
  else
    gate(k) = 1;
  end
end

r13 = nan(numel(z), 1);
r13_reason = nan(numel(z), 1);
for i = find(band).'
  k = kz(i);
  if isnan(k), continue; end
  if gate(k) ~= 0
    r13_reason(i) = gate(k);
    continue;
  end
  delta2 = 2 * gpd1 * (cum_dl(k) + dlam_int(k) * (z(i) - edges(k)));
  if cos(delta2) > r13_cos_max
    r13_reason(i) = 4;
    continue;
  end
  r13(i) = H_r13(fr.dP_hh(i, :), fr.psi, mod(-th0_int(k), pi));
  r13_reason(i) = 5 * isnan(r13(i));
end
r13_reason_key = {'resolved', 'coaxial veto', 'undefined axis', ...
  'demoted unreconciled', 'off anti-phase', 'nodes unresolvable'};

% expand to rows
th_row = nan(numel(z), 1); r_row = nan(numel(z), 1);
for k = 1:Nint
  m = kz == k;
  th_row(m) = th(k); r_row(m) = rdb(k);
end

out = struct('z_int', zc_int, 'theta_int', th, 'r_db_int', rdb, ...
  'dlam_int', dlam_int, 'theta', th_row, 'r_db', r_row, 'r13_db', r13, ...
  'theta0_int', th0_int, 'theta0_res', th0_res, ...
  'r13_reason', r13_reason, 'r13_reason_key', {r13_reason_key}, ...
  'r13_gate', gate, ...
  'J0', J0, 'J', J, 'exitflag', ex, 'n_cycles_used', n_used, ...
  'cycle_best', cyc_best, 'cycle_worsened', worsened, ...
  'edges', edges, 'fitset', fitset);

end

% ---------------------------------------------------------------- helpers
function J = H_cost(th_all, r_all, dlam_int, edges, zsub, psi, fwd, obs_n, nrm, w)
NL = numel(th_all);
layers = struct('top_m', num2cell(edges(1:NL).'), ...
  'dlam', num2cell(dlam_int.'), 'theta', num2cell(th_all.'), ...
  'r_db', num2cell(r_all.'));
mod_ = ptt.fujitaModel(layers, zsub, psi, fwd);
J = 0;
flds = {mod_.phi, mod_.dP_hh, mod_.dP_hv};
for t = 1:3
  if w(t) == 0, continue; end
  d = nrm{t}(flds{t}) - obs_n{t};
  ok = isfinite(d);
  J = J + w(t) * sum(d(ok).^2) / max(nnz(ok), 1);
end
end

function v = H_place(v, idx, p)
v(idx) = p;
end

function [mu, s] = deal_std(A)
a = A(isfinite(A));
mu = mean(a); s = std(a); if s <= 0, s = 1; end
end

function rdb = H_r13(row, psi, ax)
%H_R13 Eq. (13): r = 1/tan^2(AD/2) from the two co-pol node azimuths.
% Nodes = local minima (on the periodic row) at least 3 dB below the row
% median; exactly two are required, else NaN. That 3 dB is a NODE-DETECTION
% threshold only - it says two minima are deep enough to locate, not that
% eq. (13) is valid here. Validity is the caller's anti-phase gate.
% The pair splits the periodic azimuth into two arcs; AD is the one
% CONTAINING THE AXIS `ax` (sweep frame), because the nodes straddle the
% axis - eq. (13) on the other arc would return 1/r.
rdb = NaN;
ok = isfinite(row);
if nnz(ok) < 8, return; end
med = median(row(ok));
prv = circshift(row(:), 1).'; nxt = circshift(row(:), -1).';
ismin = row < prv & row < nxt & row < med - 3;
idx = find(ismin & ok);
if numel(idx) ~= 2, return; end
p1 = psi(idx(1)); p2 = psi(idx(2));           % p1 < p2 on [0, pi)
in_arc = ax >= p1 && ax < p2;                 % axis inside [p1, p2)?
ad = p2 - p1;
if ~in_arc, ad = pi - ad; end                 % take the arc holding the axis
ad = max(min(ad, pi - 1e-6), 1e-6);
rdb = 20 * log10(1 / tan(ad / 2)^2);
end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function v = H_const(o, fr, f, d)
%H_CONST Resolve a forward constant that ptt.ershadiFabric already burned
% into the accepted dlam. fr records what that call used, so fr wins and an
% opts value that disagrees is a hard error rather than a silent override:
% dlam is never re-fit, so running the fabric step at fc = 300 MHz and
% leaving this at the 750 MHz default would rescale every modelled phase by
% 2.5x with nothing in the cost to reveal it. opts still supplies the value
% when fr predates this recording.
has_o = isstruct(o) && isfield(o, f) && ~isempty(o.(f));
has_f = isstruct(fr) && isfield(fr, f) && ~isempty(fr.(f));
if has_o && has_f && abs(o.(f) - fr.(f)) > 1e-9 * max(abs(fr.(f)), 1)
  error('ptt:ershadiInverse:constMismatch', ...
    ['opts.%s = %g disagrees with the %g that ptt.ershadiFabric used; ' ...
     'the accepted dlam carries that constant, so the two calls must ' ...
     'agree - drop the opts override or re-run the fabric step'], ...
    f, o.(f), fr.(f));
end
if has_f, v = fr.(f); elseif has_o, v = o.(f); else, v = d; end
end
