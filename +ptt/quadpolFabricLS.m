function out = quadpolFabricLS(S, z, opts)
%QUADPOLFABRICLS Fabric from a local least-squares fit of the coherence field.
%
% out = ptt.quadpolFabricLS(S, z, opts)
%
% Fits the measured complex HH-VV coherence C(psi, z), over ALL synthetic
% azimuths at once, to the exact single-column birefringence model in a
% sliding depth window. Three physical parameters per window - the axis
% azimuth theta0, the accumulated two-way phase delta0 at the window
% centre, and its local rate ddelta/dz - plus linear nuisances: a scale
% gamma for the coherence floor, and a three-component antenna-frame
% PEDESTAL (see below). dlam is ddelta/dz converted by the same constant
% ptt.quadpolFabric uses.
%
% WHY THIS EXISTS. ptt.ershadiFabric follows the published chain and takes
% the fabric axis from the cross-polarized power minimum. On this system
% the cross-polarized channels sit on a flat instrument pedestal (cross/co
% ~ -3.6 dB with 0.35 dB of depth ripple, where fabric cross-pol must null
% through every delta = 2*pi*m), so that minimum is antenna-locked: 89.6
% +- 1.9 deg across 1790 blocks spanning all headings. This estimator
% takes theta0 from the SHAPE of the coherence field instead, which the
% co-polarized channels dominate, and models rather than gates the
% coherence structure the fabric itself creates.
%
% THE MODEL. For a column with principal axes at theta0 and accumulated
% two-way eigenmode phase difference delta(z), the deterministic coherence
% at synthetic azimuth psi is, with mu = cos 2(psi - theta0),
%
%   num = (1-mu^2)/2 + ((1+mu^2)/2) cos delta + i mu sin delta
%   den = 1 - ((1-mu^2)/2) (1 - cos delta)        (= |hh|^2 = |vv|^2)
%   C_model = gamma * num / den + pedestal(psi)
%
% within a window delta(z) = delta0 + ddelta*(z - zc).
%
% THE PEDESTAL. A reciprocal antenna-fixed cross-pol term L carried
% through the synthesis adds 2cs(<L F*> - <F L*>) to <T_hh T_vv*> - an
% odd, complex sin(2psi) term - and -4c^2s^2 <|L|^2>, an even real
% sin^2(2psi) dip. Measured directly on frame 009: at every even-2pi
% crossing, where the fabric's own signature vanishes, the coherence
% field collapses to exactly this shape (DC 0.09, 4psi harmonic 0.16-0.20)
% while the on-axis phase gradient stays clean - and the |H| = 1 fabric
% model, unable to represent it, doubled ddelta there ("regular jumps").
% So the model carries pedestal(psi) = (a1 + i a2) sin 2psi + a3 sin^2
% 2psi; all three enter linearly and are solved jointly with gamma in
% closed form at every grid node.
%
% WHY THE PEDESTAL IS ESTIMATED PER FRAME, NOT PER WINDOW. Within one
% window the pedestal coefficients are depth-constant - and so is the
% leading fabric phasor when the window sits at small ddelta*L, so for a
% fabric axis near 45 deg to the antenna axes (where mu = sin 2psi and
% the two azimuthal shapes coincide) the two are confounded per window:
% measured 2.4x inflation of the window dlam error at d = 45 vs d = 15 on
% synthetics. Across the COLUMN they are orthogonal at every theta0: the
% pedestal is an instrument constant while the fabric coefficients sweep
% sin/cos of delta0 through many fringes, averaging to zero. Default mode
% 'frame' therefore runs two passes: pedestal free per window, then the
% robust frame-level median subtracted and each window refit with only a
% strongly-ridged DEVIATION from it (anchored, not pinned - see
% H_fit_windows for why hard pinning was worse). The fitted frame
% pedestal is returned in out.pedestal: a per-frame CALIBRATION
% MEASUREMENT of the same term that antenna-locks the published chain.
%
% AMBIGUITY AND SIGN. (theta0 + 90, -delta0, -ddelta) produces the same
% field, so ddelta >= 0 is enforced and theta0 is unique modulo 180: it is
% the axis along which the (deramped) phase difference grows with depth,
% the same convention ershadiFabric reaches through the polarity of Psi.
% An isotropic window has ddelta = 0 INSIDE the parameter space, so noise
% averages to zero instead of folding to a positive floor the way |Psi|
% does.
%
% TWO-PASS USE ACROSS BLOCKS. theta0 and the pedestal vary slowly (or not
% at all) along track, so estimate both once per frame, then pass them
% back in for the per-block section (opts.theta0 and opts.pedestal set):
% the block fits then have two free physical parameters, no
% block-to-block axis jitter, and a consistent calibration.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels hh, vv, hv, vh; OR a
%         struct with field M ([Nt x 4 x 4] moments from
%         ptt.quadpolMoments / ptt.rotateMoments) when the caller has
%         already built - and possibly heading-compensated - the moment
%         matrix, as the curved-line path does
%   z     [Nt x 1] depth [m]
%   opts  fc (750e6), psi_step_deg (2), psi_offset_deg (0),
%         win_short_m (10, coherence
%         multilook), win_fit_m (60, LS window), step_m (15, window
%         centres), dlam_max (0.25), deramped (true),
%         theta0 ([] = estimate; else radians, scalar or [Nw x 1] or a
%         struct('z', 'theta') interpolated as a DOUBLED-ANGLE PHASOR),
%         pedestal ('frame' default: two-pass, frame-level estimate;
%         'window': free per window; a [1 x 3] vector [a1 a2 a3] to
%         anchor, e.g. the frame estimate handed to per-block fits; or a
%         complex [Npsi x 1] pedestal FIELD to subtract - the curved-line
%         path builds that field from the survey-calibrated antenna-frame
%         pedestal mixed over the measured heading distribution, since
%         rotating a sub-block's moments by -h maps the antenna bases
%         sin2(phi) to sin2(psi - h)),
%         q_min (0.05), dlam_min_theta (0.01), theta_step_deg (3),
%         weighting ('crb' default | 'uniform'),
%         theta_grid ([] = the full 0..180 grid at theta_step_deg; else
%         an explicit grid in radians, a vector for all windows or
%         [Nw x Ng] per window - resampling replicates search +-15 deg
%         around the full-data axis; q_theta measured over a narrow span
%         is not on the full-range scale, so such a pass must be given
%         window_ok and gates nothing itself),
%         window_ok ([] = decide per window from q_min / dlam_min_theta as
%         usual; else an [Nw x 1] logical from the full-data, full-grid
%         fit - out.window_ok - naming the windows that own an axis. A
%         resampling replicate inherits that verdict instead of
%         re-deriving it from its narrow grid, so it abstains on exactly
%         the windows its parent rejected and is never dropped for the
%         width of the grid it was given),
%         theta_const (false; true holds ONE axis for the whole column,
%         voted from the per-window theta curves - see the
%         CONSTANT-ORIENTATION MODE block), theta_const_q_min (0.02, the
%         normalised curve range below which a window does not vote),
%         z_valid (Inf: the depth [m] the record is ICE to - the bed less
%         a margin. A window reaching below it is not fitted at all: it
%         returns NaN everywhere, so it neither votes on a held axis nor
%         enters the frame pedestal median nor hands an axis to anything
%         downstream. See SUB-BED WINDOWS below.)
%
% SUB-BED WINDOWS. The record usually runs past the bed, and below it the
% channels hold noise plus the antenna-fixed pedestal and nothing else. Such
% a window still has a best-fitting axis - wherever the pedestal and the
% noise put it - and nothing in its own statistics marks it as meaningless:
% noise easily clears q_min and, because the free axis may flip a quarter
% turn to keep the rate positive, dlam_min_theta too. MEASURED on the
% constant-orientation products (23 Sep 2026): of the windows that passed
% the vote gates, 95% at Eastwind (bed ~225 m) and McMurdo (~268 m) and
% 35% at Taylor Dome (~837 m) lay below the bed, and because the vote
% pools curves normalised by their own data power each counted as much as
% an ice window. The held axis landed a median 51 / 41 / 9 deg off the
% above-bed free axis (Ridge A, whose record ends in ice: 0.4 deg). Only
% the caller knows where the bed is, so it says so here.
%
% Output fields (per window centre zw [Nw x 1])
%   zw, theta0, dlam, gamma, leak (pass-A per-window pedestal magnitude,
%   a calibration diagnostic), resid (weighted rms), q_theta, delta0,
%   pedestal ([1 x 3] frame estimate, NaN when not estimated), and
%   dlam_z / theta0_z interpolated back onto z. grad_per_dlam echoes the
%   conversion constant. theta_cost [Nw x Ngrid] is every window's theta
%   cost curve from the free grid search over theta_grid (raw weighted
%   misfit, minimised over delta0 and ddelta at each node) with theta_c2
%   its data-power normaliser - the curves the constant-orientation vote
%   pools, exported so a pooling rule can be examined offline; NaN when
%   the axis was caller-supplied. window_ok [Nw x 1] logical is the
%   per-window verdict this fit reached (free mode: which windows own an
%   axis; constant mode: which windows voted), to be handed back in as
%   opts.window_ok by a resampling pass. Constant mode adds theta_const (the held
%   axis, NaN if the vote failed and the column fell back to per-window
%   axes), theta_const_q (mean per-window contrast at that axis),
%   theta_const_n (voting windows) and theta_spread_deg (the assumption
%   check - see the block).

if nargin < 3, opts = struct(); end
fc = H_opt(opts, 'fc', 750e6);
psi_step = H_opt(opts, 'psi_step_deg', 2);
psi_off = H_opt(opts, 'psi_offset_deg', 0);
win_short = H_opt(opts, 'win_short_m', 10);
win_fit = H_opt(opts, 'win_fit_m', 60);
step_m = H_opt(opts, 'step_m', 15);
dlam_max = H_opt(opts, 'dlam_max', 0.25);
deramped = H_opt(opts, 'deramped', true);
theta0_in = H_opt(opts, 'theta0', []);
ped_in = H_opt(opts, 'pedestal', 'frame');
q_min = H_opt(opts, 'q_min', 0.05);
dlam_min_theta = H_opt(opts, 'dlam_min_theta', 0.01);
% Whether a window carries enough azimuthal contrast to own an axis is a
% property of the WINDOW, settled once by the full-data fit on the full
% grid. A resampling replicate is handed that decision instead of
% re-deriving it from the narrow grid it was given.
window_ok_in = H_opt(opts, 'window_ok', []);
% Constant fabric orientation over the fitted column. Off by default so
% existing results are unchanged; see the CONSTANT-ORIENTATION MODE block.
theta_const = H_opt(opts, 'theta_const', false);
% a window votes on the pooled axis only if its own theta cost curve varies
% by more than this fraction of its data power - a flat curve has no axis
pool_q_min = H_opt(opts, 'theta_const_q_min', 0.02);
th_step = H_opt(opts, 'theta_step_deg', 3);
% An explicit theta search grid (radians; a vector shared by every window,
% or [Nw x Ng] per window) replaces the full 0..180 grid. Resampling
% replicates use it to search only +-15 deg around the full-data axis,
% which is what makes a jackknife of the frame pass affordable.
th_grid_in = H_opt(opts, 'theta_grid', []);
% the depth the record is ice to; windows reaching below it are not fitted
z_valid = H_opt(opts, 'z_valid', Inf);

z = z(:);
Nt = numel(z);
dz = median(abs(diff(z)));

C = ptt.constants();
n_ice = sqrt(C.eps_bar);
% Same constant as ptt.quadpolFabric, cross-checked there against the
% independent fringe_check.m calibration to 0.7%.
grad_per_dlam = 2 * pi * fc * C.deps / (n_ice * C.c * 1e9);   % rad/m

% --- multilooked coherence over the full sweep. quadpolMoments averages
% the traces; the range term supplies the short window, so A.chhvv is
% already the windowed coherence of Ershadi eq. (7) at every azimuth.
if isfield(S, 'M') && ~isfield(S, 'hh')
  % caller-supplied (possibly geographic-frame) moments; the short-window
  % multilook then happens on the moment rows
  M = S.M;
  nr = max(3, round(win_short / max(dz, eps)));
  kr = ones(nr, 1) / nr;
  for k = 1:4
    for l = 1:4
      M(:, k, l) = conv(M(:, k, l), kr, 'same');
    end
  end
else
  nr = max(3, round(win_short / max(dz, eps)));
  M = ptt.quadpolMoments(S, [nr size(S.hh, 2)]);
end
% psi_offset_deg shifts the whole synthesis grid; two runs at half
% density and complementary offsets give interleaved azimuth halves for
% split-sample systematics tests.
psi = (psi_off:psi_step:psi_off+180-psi_step) * pi/180;
A = ptt.quadpolAzimuth(M, psi);
Cm = A.chhvv;
if deramped
  Cm = conj(Cm);
end

% Inverse-variance weights from the measured coherence itself: the
% Cramer-Rao variance of a coherence estimate goes as (1-|C|^2)/(2N|C|^2),
% so each (psi, z) sample carries |C|^2/(1-|C|^2), capped so one pristine
% sample cannot own a window. This suppresses the incoherent part of the
% crossing-depth corruption; the coherent part is what the pedestal terms
% are for.
if strcmpi(H_opt(opts, 'weighting', 'crb'), 'crb')
  Wc = abs(Cm).^2 ./ max(1 - abs(Cm).^2, 0.02);
  Wc = min(Wc, 25);
else
  Wc = ones(size(Cm));
end

% --- window centres and in-window sample decimation. ~2 m sampling keeps
% the fit over-determined without dragging 140 correlated samples through
% every grid evaluation.
half = win_fit / 2;
zw = (z(1) + half:step_m:z(end) - half).';
Nw = numel(zw);
jdec = max(1, round(2 / max(dz, eps)));

% theta0 supplied by the caller, in whichever of the three forms
if ~isempty(window_ok_in)
  window_ok_in = logical(window_ok_in(:));
  if numel(window_ok_in) ~= Nw
    error('ptt:quadpolFabricLS:windowOk', ...
      'opts.window_ok has %d entries for %d windows', numel(window_ok_in), Nw);
  end
end

th_fix = nan(Nw, 1);
if ~isempty(theta0_in)
  if isstruct(theta0_in)
    % Interpolate the DOUBLED-ANGLE PHASOR, never an unwrapped angle: an
    % unwrap over gappy, noisy axis samples can slip a branch, which is a
    % silent 90 deg axis error handed to every block below the slip - the
    % 0.25-railed windows at 957 m on frame 007 were exactly that.
    % Windows outside the profile's z-range take the end values: a linear
    % extrapolation of a phasor can pass near zero, and the arbitrary
    % angle it lands on would be PINNED by the block fit.
    zt = theta0_in.z(:);
    ph = interp1(zt, exp(2i*theta0_in.theta(:)), ...
      min(max(zw, min(zt)), max(zt)), 'linear');
    th_fix = 0.5 * angle(ph);
  elseif isscalar(theta0_in)
    th_fix(:) = theta0_in;
  else
    th_fix = theta0_in(:);
    if numel(th_fix) ~= Nw
      error('ptt:quadpolFabricLS:theta0', ...
        'opts.theta0 has %d entries for %d windows', numel(th_fix), Nw);
    end
  end
end

% --- everything the window fitter needs, bundled once.
% A caller-supplied grid may be a RESTRICTED search (the jackknife's
% +-15 deg around the full-data axis). q_theta measured across it is NOT
% comparable to one measured across the full 0-180 range, and nothing here
% tries to make it so: a narrow-grid pass does not adjudicate whether its
% windows have an axis at all. It is handed that verdict in
% opts.window_ok, so the contrast it computes is a diagnostic only and the
% grid's width cannot gate anything.
if isempty(th_grid_in)
  th_grid = (0:th_step:180-th_step) * pi/180;
elseif isvector(th_grid_in)
  th_grid = th_grid_in(:).';
else
  th_grid = th_grid_in;
  if size(th_grid, 1) ~= Nw
    error('ptt:quadpolFabricLS:thetaGrid', ...
      'opts.theta_grid has %d rows for %d windows', size(th_grid, 1), Nw);
  end
end
P = struct('z', z, 'zw', zw, 'half', half, 'jdec', jdec, 'psi', psi, ...
  'th_grid', th_grid, ...
  'dd_grid', linspace(0, dlam_max * grad_per_dlam, 26), ...
  'd0_grid', (0:15:345) * pi/180, ...
  'dd_max', dlam_max * grad_per_dlam, 'th_fix', th_fix, 'z_valid', z_valid);

pedestal = nan(1, 3);
if isnumeric(ped_in)
  % A caller-supplied pedestal must match one of the two documented
  % shapes exactly; anything else is a bug at the call site (most likely
  % a psi-grid mismatch with a precomputed field) and silently falling
  % back to the two-pass frame estimate would hide it.
  if numel(ped_in) == 3
    pedestal = ped_in(:).';
    R = H_fit_windows(Cm, Wc, P, pedestal);
    leak_diag = repmat(norm(pedestal), Nw, 1);
  elseif numel(ped_in) == numel(psi)
    % precomputed pedestal field over the sweep azimuths
    R = H_fit_windows(Cm, Wc, P, ped_in(:));
    leak_diag = repmat(max(abs(ped_in(:))), Nw, 1);
  else
    error('ptt:quadpolFabricLS:pedestal', ...
      ['numeric opts.pedestal must be [1 x 3] coefficients or a field ' ...
      'over the %d sweep azimuths; got %d entries'], ...
      numel(psi), numel(ped_in));
  end
elseif (ischar(ped_in) || isstring(ped_in)) && strcmpi(ped_in, 'window')
  R = H_fit_windows(Cm, Wc, P, []);
  leak_diag = R.leak;
  okp = all(isfinite(R.ped_coef), 2);
  if nnz(okp) >= 3
    pedestal = median(R.ped_coef(okp, :), 1);
  end
elseif (ischar(ped_in) || isstring(ped_in)) && strcmpi(ped_in, 'frame')
  % the default two-pass
  Ra = H_fit_windows(Cm, Wc, P, []);
  leak_diag = Ra.leak;
  okp = all(isfinite(Ra.ped_coef), 2);
  if nnz(okp) >= 5
    pedestal = median(Ra.ped_coef(okp, :), 1);
    R = H_fit_windows(Cm, Wc, P, pedestal);
  else
    R = Ra;   % too few windows to trust a frame estimate
  end
else
  error('ptt:quadpolFabricLS:pedestal', ...
    'opts.pedestal must be ''frame'', ''window'', [1 x 3], or a field');
end

% --- CONSTANT-ORIENTATION MODE. One theta0 for the whole column instead of
% one per window: the windows vote by their own theta cost curves, which the
% grid search already computed, and the winner is refit into every window.
%
% WHY IT IS OFTEN THE BETTER MODEL. Per-window theta0 spends one free
% parameter per window on a quantity that at a divide does not vary: Ridge
% A's frame axis is 96.5 +- 7.4 deg across 36 frames, while its per-window
% scatter within a frame is ~13 deg - i.e. most of that scatter is
% estimation noise, not ice. Removing those degrees of freedom leaves the
% depth-varying quantity that IS physical (delta0 and its gradient, hence
% dlam) better determined, and hands downstream steps an axis stable enough
% to condition on.
%
% WEIGHTING. Each window's curve is NORMALISED by its own data power
% (theta_c2) before the sum, i.e. pooled on the q_theta scale, so a window
% votes in proportion to how sharply ITS curve is peaked and nothing else.
% Raw-cost pooling was tried first and is wrong on real data: the data
% power - the CRB-weighted coherence power - runs 100-1000x higher in the
% near-surface windows than at depth, and on Ridge A frame 009 seven
% windows at 57-147 m carried 86% of the raw vote and put the axis at 102
% deg while every window from 200 m to the bottom sat at 93-96 deg (the
% frame's published value). Those surface windows are exactly where the
% co-pol reference offset makes the model untrustworthy, and a raw sum
% lets them outvote the whole column. Normalised, the vote lands on 96
% with or without any depth cut, so no site-specific depth enters here.
% A flat curve still carries nothing - its normalised range is small - and
% windows flat within pool_q_min of their own power are excluded outright.
%
% WHEN NOT TO USE IT. It is a MODEL ASSUMPTION, not a refinement: where the
% axis genuinely rotates with depth - Thwaites' margin, EastGRIP - it will
% return some average of the rotation and the per-window residuals will say
% so. Check out.theta_const_q (the pooled contrast) and the residual
% profile before adopting it at a new site.
th_c = NaN; theta_const_q = NaN; theta_spread_deg = NaN; n_theta_pool = 0;
held = false; vote_ok = [];
Rfree = R;   % the per-window pass: its curves ARE the vote, and stay the diagnostic
if theta_const && isempty(theta0_in)
  if size(P.th_grid, 1) > 1
    error('ptt:quadpolFabricLS:thetaGrid', ...
      'constant-theta mode pools curves across windows and needs one shared theta_grid');
  end
  ok_w = all(isfinite(R.cost_th), 2) & isfinite(R.c2) & R.c2 > 0;
  if ~isempty(window_ok_in)
    % a replicate pools exactly the windows its parent accepted, so every
    % replicate votes on the same set and the pooled axis is comparable
    % across them. Re-deciding here from a narrow grid is what left
    % replicates voting on different window sets, and a segment where
    % enough of them fell below the floor returned theta_const = NaN.
    ok_w = ok_w & window_ok_in;
  else
    % a window only votes if its own curve is not flat ...
    rng_w = max(R.cost_th, [], 2) - min(R.cost_th, [], 2);
    ok_w = ok_w & rng_w > pool_q_min * max(R.c2, realmin);
    % ... and only if the free fit would have TRUSTED its axis: the same
    % q_min / dlam_min_theta rule that abstains a per-window theta0 below.
    % A window with no phase gradient has a minimum, but it is where noise
    % put it - on frame 009 the 160-190 m windows, where the co-pol power
    % collapses, all bottomed a quarter turn off (the (theta+90, -delta)
    % ambiguity) with q ~0.5. Let in, they do little to the vote but they
    % turn the spread diagnostic into a false rotation.
    ok_w = ok_w & R.q_theta >= q_min & R.dlam / grad_per_dlam >= dlam_min_theta;
  end
  vote_ok = ok_w;
  if nnz(ok_w) >= 3
    Cn = R.cost_th(ok_w, :) ./ max(R.c2(ok_w), realmin);   % q_theta scale
    pooled = sum(Cn, 1);
    [~, ib] = min(pooled);
    th_c = P.th_grid(ib);
    % Sub-grid refinement: a parabola through the minimum and its two
    % neighbours. Without it the held axis is quantised to the grid step
    % (3 deg by default) and so is anything derived from it - a jackknife
    % over sub-blocks reports zero spread whenever every replicate lands
    % on the same node, which is not an uncertainty of 0 deg. The full
    % grid is periodic in pi, so its end nodes wrap; a reduced grid
    % (theta_grid supplied) has real edges and is left alone there.
    Ng = numel(pooled);
    full_grid = isempty(th_grid_in);
    if Ng >= 3 && (full_grid || (ib > 1 && ib < Ng))
      im = mod(ib - 2, Ng) + 1; ip = mod(ib, Ng) + 1;
      pm = pooled(im); p0 = pooled(ib); pp = pooled(ip);
      curv = pm - 2*p0 + pp;
      if curv > 0
        step_g = abs(angle(exp(2i * (P.th_grid(min(ib+1, Ng)) - P.th_grid(max(ib-1, 1)))))) / (2 * (min(ib+1, Ng) - max(ib-1, 1)));
        th_c = th_c + 0.5 * (pm - pp) / curv * step_g;
      end
    end
    % pooled contrast: the mean per-window contrast at the pooled axis, on
    % the same scale q_theta uses per window
    theta_const_q = (max(pooled) - min(pooled)) / nnz(ok_w);
    % IS THE ASSUMPTION TRUE HERE? The pooled contrast does NOT answer that
    % - measured on synthetics it came out HIGHER on a rotating column
    % (0.225) than on a constant one (0.219), because it reports how
    % sharply the pooled curve is peaked, not whether the windows agreed.
    % What answers it is the SPREAD of the per-window minima: windows of a
    % constant column all bottom at the same theta, windows of a rotating
    % one bottom across the rotation. Reported as a circular standard
    % deviation of the doubled angle, in degrees, with each minimum
    % weighted by its window's contrast: a near-flat window's minimum is
    % where noise put it (on frame 009 four such windows sat a quarter
    % turn off, the (theta+90, -delta) ambiguity), and unweighted they
    % would report a rotation the column does not have.
    [~, iw] = min(R.cost_th(ok_w, :), [], 2);
    th_w = P.th_grid(iw);
    qw = max(R.q_theta(ok_w), 0);
    Rw = sum(qw(:) .* exp(2i * th_w(:))) / max(sum(qw), realmin);
    theta_spread_deg = rad2deg(sqrt(max(-2*log(max(abs(Rw), realmin)), 0))) / 2;
    P.th_fix = repmat(th_c, Nw, 1);
    if isnan(pedestal(1))
      R = H_fit_windows(Cm, Wc, P, []);
    else
      R = H_fit_windows(Cm, Wc, P, pedestal);
    end
    R.theta0(:) = th_c;                 % held, not re-estimated
    % ... except where no window was fitted at all - below z_valid, or
    % with too little finite data - which is not a weak window the column
    % speaks for but no measurement, and must not read as one downstream
    R.theta0(~isfinite(R.c2)) = NaN;
    % The held refit computes no theta curve (nothing to search), so the
    % per-window axis information and the vote it was drawn from come
    % from the free pass: q_theta stays "how much axis does THIS window
    % carry", which the handoff and the spread diagnostic both need.
    R.q_theta = Rfree.q_theta;
    R.cost_th = Rfree.cost_th;
    R.c2 = Rfree.c2;
    n_theta_pool = nnz(ok_w);
    held = true;
  else
    n_theta_pool = nnz(ok_w);
    warning('ptt:quadpolFabricLS:noPool', ...
      ['constant-theta mode: only %d window(s) carry a trusted, non-flat ' ...
      'theta curve, need 3; falling back to per-window theta0'], nnz(ok_w));
  end
end

theta0 = R.theta0; dlam = R.dlam / grad_per_dlam; gam = R.gam;
resid = R.resid; q_theta = R.q_theta; delta0 = R.delta0;

% Per-window abstention on the axis: a window with no fabric or a flat
% curve has no axis of its own. A HELD column is exempt for the same
% reason a caller-supplied theta0 is - the axis is the column's, not the
% window's, and dropping it from the near-isotropic windows would hand
% those depths back to a different model downstream (the blocks would
% fall back to the frame profile there, mixing held and free axes in one
% section).
% A resampling replicate does not re-decide this. It is handed the parent
% window's verdict in opts.window_ok and abstains on exactly the windows
% the parent rejected, whatever grid it searched - its own contrast, which
% is not on the full-range scale, gates nothing. Re-deriving it from a
% +-15 deg grid is what dropped replicate values for having been
% narrow-searched, silently, since nothing downstream reports it.
if isempty(theta0_in) && ~held
  if isempty(window_ok_in)
    window_ok = q_theta >= q_min & dlam >= dlam_min_theta;
  else
    window_ok = window_ok_in;
  end
  theta0(~window_ok) = NaN;
elseif ~isempty(vote_ok)
  % held: the axis is the column's, so the windows that DECIDED it are the
  % set a replicate must inherit, not the ones that carry a theta0
  window_ok = vote_ok;
else
  window_ok = isfinite(theta0);
end

% interpolate back onto z for the section plumbing; the phasor keeps the
% axis average honest across the modulo-180 wrap
ok = isfinite(theta0);
theta0_z = nan(Nt, 1); dlam_z = nan(Nt, 1);
if nnz(ok) >= 2
  ph = interp1(zw(ok), exp(2i*theta0(ok)), z, 'linear');
  theta0_z = mod(angle(ph)/2, pi);
end
% dlam keeps its abstention NaNs: interpolating over the full window grid
% fills between measured neighbours but never bridges across an abstained
% window, so the section (and the figure, which paints NaN white) can tell
% abstained cells from measured ones
okd = isfinite(dlam);
if nnz(okd) >= 2
  dlam_z = interp1(zw, dlam, z, 'linear');
end

out = struct('zw', zw, 'theta0', theta0, 'dlam', dlam, 'gamma', gam, ...
  'leak', leak_diag, 'pedestal', pedestal, ...
  'resid', resid, 'q_theta', q_theta, 'delta0', delta0, ...
  'theta0_z', theta0_z, 'dlam_z', dlam_z, ...
  'theta_const', th_c, 'theta_const_q', theta_const_q, ...
  'theta_spread_deg', theta_spread_deg, ...
  'theta_const_n', n_theta_pool, ...
  'theta_grid', P.th_grid, 'theta_cost', Rfree.cost_th, ...
  'theta_c2', Rfree.c2, 'window_ok', window_ok, ...
  'grad_per_dlam', grad_per_dlam, 'psi', psi);

end

function v = H_opt(o, f, d)
if isstruct(o) && isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function R = H_fit_windows(Cm, Wc, P, ped)
%H_FIT_WINDOWS One pass of per-window fits.
%
% ped = []       pedestal free: solved with gamma at every node (4x4),
%                per-window coefficients returned in R.ped_coef
% ped = [1 x 3]  pedestal ANCHORED: the frame estimate is subtracted from
%                the data and each window fits only a DEVIATION from it,
%                under a much stronger ridge. Hard pinning was tried and
%                made things worse: the pedestal is constant in raw units
%                but its contribution to the NORMALIZED coherence varies
%                with depth and azimuth through the co-pol power in the
%                normalization, so windows at the crossings need local
%                headroom. The anchor gives it to them while still
%                breaking the 45 deg per-window degeneracy: where the
%                geometry cannot separate pedestal from fabric, the
%                deviation shrinks to zero and the frame value rules.
Nw = numel(P.zw);
R = struct('theta0', nan(Nw,1), 'dlam', nan(Nw,1), 'gam', nan(Nw,1), ...
  'resid', nan(Nw,1), 'q_theta', nan(Nw,1), 'delta0', nan(Nw,1), ...
  'leak', nan(Nw,1), 'ped_coef', nan(Nw,3), ...
  'cost_th', nan(Nw, size(P.th_grid, 2)), 'c2', nan(Nw,1));
psi = P.psi;
sb = sin(2 * psi(:));
s2b = sb.^2;
free_ped = isempty(ped);
tau_rel = 0.02;          % light ridge while the pedestal is being found
if ~free_ped
  if numel(ped) == 3
    pf = (ped(1) + 1i*ped(2)) * sb + ped(3) * s2b;   % [Nk x 1]
  else
    pf = ped(:);          % precomputed field over the sweep azimuths
  end
  tau_rel = 0.15;        % strong anchor on the per-window deviation
end
os = optimset('Display', 'off', 'MaxFunEvals', 400, 'MaxIter', 400, ...
  'TolFun', 1e-6, 'TolX', 1e-6);

for w = 1:Nw
  % a window that reaches below the ice is noise and pedestal only; see
  % SUB-BED WINDOWS in the header
  if P.zw(w) + P.half > P.z_valid, continue; end
  jj = find(P.z >= P.zw(w) - P.half & P.z <= P.zw(w) + P.half);
  jj = jj(1:P.jdec:end);
  Cw = Cm(jj, :).';                 % [Nk x Nj]: azimuth rows, depth cols
  u = (P.z(jj) - P.zw(w)).';        % [1 x Nj]
  ok = isfinite(Cw);
  if nnz(ok) < 0.5 * numel(Cw), continue; end
  Cw(~ok) = 0;
  wgt = Wc(jj, :).';
  wgt(~ok) = 0;
  if ~free_ped
    Cw = Cw - pf;                   % implicit expansion over depth
    Cw(~ok) = 0;
  end
  C2 = sum(wgt .* abs(Cw).^2, 'all');
  if C2 <= 0, continue; end

  wrow = sum(wgt, 2);
  rCrow = sum(wgt .* real(Cw), 2);
  iCrow = sum(wgt .* imag(Cw), 2);
  K = struct('W2', (sb.^2).' * wrow, 'W3', (sb.^3).' * wrow, ...
    'W4', (sb.^4).' * wrow, 'd1', sb.' * rCrow, 'd2', sb.' * iCrow, ...
    'd3', s2b.' * rCrow, 'sb', sb, 's2b', s2b);
  K.tau = tau_rel * max(K.W2, realmin);

  if isfinite(P.th_fix(w))
    ths = P.th_fix(w);
  else
    ths = P.th_grid(min(w, size(P.th_grid, 1)), :);
  end
  best = struct('cost', inf, 'th', NaN, 'd0', NaN, 'dd', NaN);
  cost_th = inf(numel(ths), 1);
  for it = 1:numel(ths)
    mu = cos(2 * (psi(:) - ths(it)));
    Ak = (1 - mu.^2) / 2;
    Bk = (1 + mu.^2) / 2;
    for id = 1:numel(P.dd_grid)
      cu = cos(P.dd_grid(id) * u); su = sin(P.dd_grid(id) * u);
      for i0 = 1:numel(P.d0_grid)
        c0 = cos(P.d0_grid(i0)); s0 = sin(P.d0_grid(i0));
        cd = c0 * cu - s0 * su;
        sd = s0 * cu + c0 * su;
        num = Ak + Bk * cd + 1i * (mu * sd);
        den = max(1 - Ak * (1 - cd), 0.05);
        H = num ./ den;
        rH = real(H); iH = imag(H);
        G1 = sum(wgt .* (real(Cw) .* rH + imag(Cw) .* iH), 'all');
        G2 = sum(wgt .* (rH.^2 + iH.^2), 'all');
        rHrow = sum(wgt .* rH, 2);
        iHrow = sum(wgt .* iH, 2);
        cost = H_solve(C2, G1, G2, sb.' * rHrow, sb.' * iHrow, ...
          s2b.' * rHrow, K);
        if cost < cost_th(it), cost_th(it) = cost; end
        if cost < best.cost
          best = struct('cost', cost, 'th', ths(it), ...
            'd0', P.d0_grid(i0), 'dd', P.dd_grid(id));
        end
      end
    end
  end
  if ~isfinite(best.cost), continue; end

  % theta0 contrast: how much worse the fit gets a quarter turn away. An
  % isotropic window is flat here and its theta0 means nothing.
  if numel(ths) > 1
    R.q_theta(w) = (max(cost_th) - best.cost) / max(C2, realmin);
    % the whole theta cost curve, kept so a CONSTANT-orientation fit can
    % pool it across windows instead of re-running the grid search
    R.cost_th(w, :) = cost_th(:).';
  end
  R.c2(w) = C2;

  % polish from the best grid node. fminsearch is base MATLAB; the linear
  % nuisances stay closed-form inside the objective, and ddelta is kept
  % in range by a clamp rather than a penalty cliff. A caller-fixed
  % theta0 stays PINNED: the whole point of the two-pass use is that the
  % block fits cannot wander off the frame axis.
  fun = @(p) H_cost(p, psi, u, Cw, wgt, C2, P.dd_max, K);
  if isfinite(P.th_fix(w))
    p2 = fminsearch(@(q) fun([P.th_fix(w); q(:)]), [best.d0; best.dd], os);
    p = [P.th_fix(w); p2(:)];
  else
    p = fminsearch(fun, [best.th; best.d0; best.dd], os);
  end
  [cst, g, la, ac] = fun(p);
  dd = min(max(p(3), 0), P.dd_max);

  R.theta0(w) = mod(p(1), pi);
  R.delta0(w) = mod(p(2), 2*pi);
  if dd > 0.98 * P.dd_max
    % Railed at the cap: the window found no interior optimum, so the
    % value is the bound, not a rate. Abstain rather than report it.
    R.dlam(w) = NaN;
  else
    R.dlam(w) = dd;   % rad/m; the caller converts by grad_per_dlam
  end
  R.gam(w) = g;
  R.leak(w) = la;
  R.ped_coef(w, :) = ac;
  R.resid(w) = sqrt(max(cst, 0) / C2);
end
end

function [cost, gamma, amp, ac] = H_cost(p, psi, u, Cw, wgt, C2, dd_max, K)
% Weighted misfit of the model at (theta0, delta0, ddelta), with the
% linear amplitudes eliminated in closed form: gamma alone when the
% pedestal is pinned (K empty), gamma plus the three pedestal terms
% otherwise (H_solve).
th = p(1);
d0 = p(2);
dd = min(max(p(3), 0), dd_max);
mu = cos(2 * (psi(:) - th));
Ak = (1 - mu.^2) / 2;
Bk = (1 + mu.^2) / 2;
cd = cos(d0 + dd * u);
sd = sin(d0 + dd * u);
num = Ak + Bk * cd + 1i * (mu * sd);
den = max(1 - Ak * (1 - cd), 0.05);
H = num ./ den;
rH = real(H); iH = imag(H);
G1 = sum(wgt .* (real(Cw) .* rH + imag(Cw) .* iH), 'all');
G2 = sum(wgt .* (rH.^2 + iH.^2), 'all');
rHrow = sum(wgt .* rH, 2);
iHrow = sum(wgt .* iH, 2);
[cost, gamma, amp, ac] = H_solve(C2, G1, G2, K.sb.' * rHrow, ...
  K.sb.' * iHrow, K.s2b.' * rHrow, K);
end

function [cost, gamma, amp, ac] = H_solve(C2, G1, G2, b01, b02, b03, K)
% Joint closed-form solve for the model scale gamma and the pedestal
% amplitudes (a1 + i a2 on sin 2psi, a3 on sin^2 2psi): all four enter
% the model linearly, so the weighted LS reduces to a 4x4 normal system.
% gamma is clamped to [0, 1.05] (a coherence above 1 is not physical) by
% refitting the nuisances at the clamped value, so a clamp never buys a
% better cost than the data support. The ridge tau keeps the nuisances
% quiet where a single window cannot separate them from the fabric term;
% the frame-level median across windows is what separates them for good.
M = [G2,  b01,          b02,          b03; ...
     b01, K.W2 + K.tau, 0,            K.W3; ...
     b02, 0,            K.W2 + K.tau, 0; ...
     b03, K.W3,         0,            K.W4 + K.tau];
r = [G1; K.d1; K.d2; K.d3];
x = M \ r;
if ~all(isfinite(x))
  cost = C2; gamma = 0; amp = 0; ac = nan(1, 3);
  return;
end
if x(1) < 0 || x(1) > 1.05
  g = min(max(x(1), 0), 1.05);
  xa = M(2:4, 2:4) \ (r(2:4) - g * M(2:4, 1));
  x = [g; xa];
end
cost = C2 - 2 * (x.' * r) + x.' * M * x;
gamma = x(1);
ac = x(2:4).';
amp = norm(ac);
end
