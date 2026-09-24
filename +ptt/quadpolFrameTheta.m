function fp = quadpolFrameTheta(T, z, az_tr, x_along, opts)
%QUADPOLFRAMETHETA Frame pass of the quad-pol chain, laterally segmented.
%
% fp = ptt.quadpolFrameTheta(T, z, az_tr, x_along, opts)
%
% Estimates the frame-level pedestal and the GEOGRAPHIC fabric-axis
% profile theta0(z) that the per-block section fits inherit - the "first
% pass" of the two-pass design - and, when the frame is long enough,
% re-estimates theta0(z) per ALONG-TRACK SEGMENT so that lateral axis
% variation is resolved instead of pooled away.
%
% WHY. The two-pass design's single point of failure is the frame pass
% pooling every trace into one theta0(z) handoff: an along-frame change in
% the fabric (EastGRIP's 6.1 km frame crosses a shear margin mid-way -
% 19.9 km before the pipeline rebuilt that season's trajectory)
% decoheres the pooled fit, and every block then inherits the result.
% test_egrip_blocks.m measures what that costs at the season's real 2.83 m
% spacing, and it does NOT need a dead frame to bite: at a modest 0.010
% dlam/km the pooled axis is only ~3 deg off over the 300-1100 m band,
% yet blocks handed that profile recover dlam ~12x worse there than
% blocks handed the true axis (verdict 2). Push the ramp further and the
% pooled fit abstains rather than lies - coverage falls away, and every
% rate that still carries enough band coverage to be scored keeps its
% axis within a few degrees of truth (verdict 1) - so the damage is
% upstream of the blocks either way, and no block size addresses it.
% Blocks themselves cannot self-rescue -
% 14-126 traces is too little coherence to estimate theta0 alone; ~2 km
% is plenty. So the fix is lateral resolution of the FIRST pass, not the
% second.
%
% WHAT STAYS FRAME-LEVEL. The pedestal: it is an instrument constant
% (frame-stable a3 across every surveyed frame, sign set per season), so
% one antenna-frame estimate per frame remains correct and segments
% receive it as a fixed field, never refit. The frame-pooled profile is
% also kept - it is the health signal, the nseg == 1 result, and the
% fallback for segments with nothing usable.
%
% HOW SEGMENTS FIT. Every segment fits in the GEOGRAPHIC frame from
% pooled heading-rotated sub-block moments, exactly the machinery the
% curved-line path validated on a 90 deg arc (test_quadpol_curved.m): on
% a straight frame the rotation is a constant and the path reduces to the
% straight fit. The segment's pedestal field is the frame antenna-frame
% pedestal mixed over the SEGMENT's own heading distribution
% (rotation maps sin2(phi) -> sin2(psi - h)). Segment profiles are
% smoothed in depth exactly like the frame profile (q_theta-weighted
% doubled-angle phasor, never an unwrap), then smoothed ACROSS segments
% the same way, and windows where a segment has no usable theta fall back
% to the frame profile - so a short or dead segment degrades to current
% behavior instead of poisoning its blocks.
%
% Inputs
%   T        struct of complex [Nt x Nx] channels hh, vv, hv, vh
%   z        [Nt x 1] depth [m]
%   az_tr    [Nx] per-trace geographic track AXIS (deg, mod 180), smoothed
%            upstream (the pipeline builds it from the trajectory on the
%            doubled-angle phasor)
%   x_along  [Nx] along-track position [m] (cumulative distance; any
%            origin). Segment boundaries are cut in this coordinate so
%            "2 km" means 2 km at every site regardless of trace spacing.
%   opts     fc (750e6), dlam_max (0.25), deramped (true),
%            seg_len_m (2000; frames shorter than 2 segments return the
%            frame profile unchanged), curved_p95_deg (5),
%            nblk_rot (200 traces per heading-rotation sub-block),
%            psi_step_seg_deg (2, the segment/curved fit grid),
%            track_az ([] = circular mean of az_tr; the pipeline passes
%            its endpoint-derived value so the straight frame path is
%            bit-identical to the unsegmented pipeline),
%            min_seg_windows (5, fewer usable theta windows marks the
%            segment dead),
%            theta_const (false; true holds ONE axis per segment, constant
%            in depth - see below),
%            jackknife (false; true adds sub-block jackknife standard
%            errors per segment, see ptt.quadpolJackknife),
%            z_bed ([] = no bed known; else [Nx] per-trace bed depth [m],
%            NaN where unpicked) and bed_margin_m (20): every trace is
%            blanked below its own bed less the margin (ptt.maskBelowBed),
%            so each moment averages only traces still in ice at that
%            depth, and a window is fitted only where at least HALF the
%            traces of its pass (the frame, or the segment) are still in
%            ice at the window's bottom edge (ptt.quadpolFabricLS z_valid
%            at the median bed). The two together follow a bed that varies
%            along the line - 567-996 m inside one Taylor Dome frame, 42-
%            300 m on one Eastwind line - which no single depth cut can:
%            a 5th-percentile cut, tried first, left an Eastwind frame no
%            window at all (p5 79 m against a 224 m median). Without a bed
%            the record's coherent sub-bed returns vote on the held axis
%            and feed the pedestal median as if they were ice - see
%            SUB-BED WINDOWS in ptt.quadpolFabricLS.
%            A pass's windows at or below its own z_valid are NOT ice and
%            are never filled: th_seg stays NaN and q_seg 0 there in both
%            modes, neither borrowed from a live neighbouring segment nor
%            taken from the frame profile (a held segment's axis below
%            its bed would otherwise be its neighbour's, and a block that
%            clamps to that segment would inherit a depth-varying axis
%            the constant model never asserted). The two fallbacks apply
%            only to dead windows ABOVE the pass's z_valid.
%
% Output struct fp
%   zw        [Nw x 1] window centres
%   th_frame  [Nw x 1] SMOOTHED frame-pooled geographic profile (fallback
%             and health signal; NaN where nothing usable)
%   q_frame   [Nw x 1] its q_theta weights (0 where dead)
%   th_seg    [Nw x Nseg] smoothed geographic profile per segment; every
%             finite value is usable, NaN only where even the frame
%             profile is dead or the window is at or below the segment's
%             own z_valid (never filled, see opts.z_bed)
%   q_seg     [Nw x Nseg] weights (0 where the segment fell back or the
%             window is at or below its z_valid)
%   dlam_seg  [Nw x Nseg] per-segment dlam (diagnostic; blocks still fit
%             their own)
%   resid_seg [Nw x Nseg]
%   seg_x     [1 x Nseg] segment centres in x_along coordinates
%   seg_n     [1 x Nseg] traces per segment
%   nseg, curved, hspread, track_az
%   ped_ant   [1 x 3] antenna-frame pedestal. ptt.pedestalMarker() -
%             exactly [0 0 0] - MARKS a frame whose pedestal fit did not
%             converge, whose segments were then fitted without one; it is
%             an audit marker to drop, never a measurement (a real fit
%             never lands on exact zeros). Consumers test it with
%             ptt.pedestalFailed rather than for finiteness, which the
%             marker passes.
%   lsq       the frame-pass LS output (report/save compatibility)
%   th_spread_seg   [1 x Nseg] circular spread (deg) of the per-window
%             axis minima inside each segment, constant mode only (NaN
%             otherwise): the "is the axis really constant here" check.
%             For nseg == 1 it is the frame's own value, th_spread_frame.
%   held_seg  [1 x Nseg] true where the segment's axis was held (constant
%             mode and the vote succeeded); a false entry in constant mode
%             marks a segment that fell back to per-window axes
%   se_theta_seg [Nw x Nseg], se_dlam_seg [Nw x Nseg]  jackknife standard
%             errors (radians, dlam units) of the segment's raw profile,
%             NaN unless opts.jackknife; jack_n / jack_edge [1 x Nseg] the
%             replicate count and the count that hit the theta search edge
%   se_theta_n / se_theta_r / se_theta_sat [Nw x Nseg]  what a NaN in
%             se_theta_seg means: the replicates that contributed, their
%             resultant length, and whether they scattered past what the
%             bounded axis statistic can resolve (as opposed to there
%             being too few). se_theta_seg must not be read bare - see
%             ptt.circAxisSE
%
% The per-block handoff belongs to ptt.thetaProfileAt(fp, x), which
% interpolates th_seg across segments on the doubled-angle phasor and
% returns the struct('z','theta') form ptt.quadpolFabricLS accepts.
%
% See also ptt.quadpolFabricLS, ptt.thetaProfileAt, ptt.rotateMoments.

if nargin < 5, opts = struct(); end
FC = H_opt(opts, 'fc', 750e6);
DLAM_MAX = H_opt(opts, 'dlam_max', 0.25);
DERAMPED = H_opt(opts, 'deramped', true);
SEG_LEN = H_opt(opts, 'seg_len_m', 2000);
CURV_P95 = H_opt(opts, 'curved_p95_deg', 5);
NBLK_ROT = H_opt(opts, 'nblk_rot', 200);
PSI_STEP_SEG = H_opt(opts, 'psi_step_seg_deg', 2);
track_az = H_opt(opts, 'track_az', []);
MIN_SEG_W = H_opt(opts, 'min_seg_windows', 5);
% CONSTANT AXIS WITH DEPTH, PER SEGMENT. The segment is already the lateral
% unit, so this makes the model "one orientation per place, free to change
% along the survey" - which is what a divide fabric actually does, and what
% the depth smoothing below has been approximating with a kernel instead of
% stating. Off by default: switching it on changes every fitted theta0, so
% it is a decision per survey, not a default. At a site whose axis rotates
% with depth (Thwaites' margin, EastGRIP) it is the WRONG model and the
% per-window profile should be kept - check the pooled contrast
% out.theta_const_q and the residual profile before adopting it.
TH_CONST = H_opt(opts, 'theta_const', false);
% STANDARD ERRORS by delete-one jackknife over the heading sub-blocks
% (ptt.quadpolJackknife), per segment - and for the frame when it is the
% one segment. Off by default because it costs ~n_sub reduced-grid refits
% per segment; the pipeline turns it on.
JACK = H_opt(opts, 'jackknife', false);
% grid/window pass-throughs, so tests can trade resolution for speed while
% the pipeline keeps the estimator defaults
PASS_THRU = {'psi_step_deg', 'theta_step_deg', 'win_short_m', ...
  'win_fit_m', 'step_m', 'weighting'};
base = struct('fc', FC, 'deramped', DERAMPED, 'dlam_max', DLAM_MAX);
for k = 1:numel(PASS_THRU)
  if isfield(opts, PASS_THRU{k}) && ~isempty(opts.(PASS_THRU{k}))
    base.(PASS_THRU{k}) = opts.(PASS_THRU{k});
  end
end
segbase = base;
segbase.psi_step_deg = PSI_STEP_SEG;
% The constant-axis option goes on every pass whose AXIS is used - the
% frame profile (the nseg == 1 handoff and every segment's fallback) and
% the segments - and on nothing else: the curved-frame antenna pass below
% exists for its pedestal alone, and a pooled refit there would cost time
% for a theta nobody reads.
fbase = base;
fbase.theta_const = TH_CONST;
segbase.theta_const = TH_CONST;

CHAN = {'hh', 'vv', 'hv', 'vh'};
Nx = size(T.hh, 2);
az_tr = az_tr(:).';
x_along = x_along(:).';

% --- where the ice ends. Every frame-level pass stops at the frame's bed,
% and each segment below at its own.
Z_BED = H_opt(opts, 'z_bed', []);
BED_MARGIN = H_opt(opts, 'bed_margin_m', 20);
if ~isempty(Z_BED)
  Z_BED = Z_BED(:).';
  if numel(Z_BED) ~= Nx
    error('ptt:quadpolFrameTheta:zBed', ...
      'opts.z_bed has %d entries for %d traces', numel(Z_BED), Nx);
  end
  % idempotent: a caller that already masked (the pipeline does, for its
  % section blocks) loses nothing by this second pass
  T = ptt.maskBelowBed(T, z, Z_BED, BED_MARGIN);
end
zv_frame = H_zvalid(Z_BED, 1:Nx, BED_MARGIN, Inf);
base.z_valid = zv_frame;
fbase.z_valid = zv_frame;
segbase.z_valid = zv_frame;

% --- heading statistics on the doubled-angle phasor
az0 = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr)))))/2, 180);
if isempty(track_az), track_az = az0; end
hdev = abs(mod(az_tr - az0 + 90, 180) - 90);
hspread = prctile(hdev, 95);
curved = hspread > CURV_P95;

% --- frame pass: pedestal (antenna frame, instrument constant) and the
% frame-pooled geographic profile. Identical to the unsegmented pipeline.
PSI_FIT = (0:PSI_STEP_SEG:180-PSI_STEP_SEG) * pi/180;
if ~curved
  lsq = ptt.quadpolFabricLS(T, z, fbase);
  ped_ant = lsq.pedestal;
  % The frame-mode pedestal is NaN whenever fewer than five windows gave a
  % finite coefficient - and quadpolFabricLS still returns a populated
  % theta0 in that case, so the frame pass looks healthy. Left unguarded a
  % NaN pedestal would enter H_ped_field, poison every segment's coherence
  % field, kill every segment fit on the MIN_SEG_W test, and drop the whole
  % pass back onto the frame profile: nseg > 1 and ls_theta_seg populated,
  % yet bit-identical to the pooled handoff. That silent no-op would hit
  % exactly the degraded frames segmentation exists for, so it is guarded
  % here as it already is on the curved branch. The substituted EXACT zeros
  % are the audit signal - a real fit never lands on them - and the warning
  % puts it in the frame log.
  if ptt.pedestalFailed(ped_ant)
    warning('ptt:quadpolFrameTheta:pedestalFailed', ...
      ['frame pedestal did not converge; segments fit with a zero ' ...
      'pedestal (saved ls_pedestal is exactly [0 0 0] to mark it)']);
    ped_ant = ptt.pedestalMarker();
  end
  th_geo_raw = lsq.theta0 + deg2rad(track_az);
else
  % antenna-frame single pass: theta/dlam are junk on a curve, but the
  % pedestal estimate is CLEANER than on a straight line (fabric smears,
  % instrument adds coherently)
  oa = base; oa.pedestal = 'window';
  lsa = ptt.quadpolFabricLS(T, z, oa);
  ped_ant = lsa.pedestal;
  if ptt.pedestalFailed(ped_ant), ped_ant = ptt.pedestalMarker(); end
  Mg = H_geo_moments(T, az_tr, 1:Nx, NBLK_ROT, CHAN);
  og = segbase;
  og.pedestal = H_ped_field(ped_ant, az_tr, 1:Nx, PSI_FIT);
  lsq = ptt.quadpolFabricLS(struct('M', Mg), z, og);
  th_geo_raw = lsq.theta0;   % the geographic fit reports geographically
end
zw = lsq.zw;
Nw = numel(zw);
half = lsq.win_half;
% A HELD column needs no depth smoothing - the profile is one number - and
% must not get any: the kernel's weight threshold would drop windows from
% a profile that is constant by construction. Its weight is the pooled
% contrast, uniform over the windows that were fitted, so the handoff
% interpolation between segments stays depth-independent too (a
% depth-varying weight would turn two constant neighbours into a
% depth-varying blend); a window the column never fitted (below z_valid)
% carries no axis and so no weight. A column whose vote failed
% (theta_const NaN) fell back to per-window axes inside the estimator and
% is treated as such here.
th_spread_frame = NaN;
if isfield(lsq, 'theta_spread_deg'), th_spread_frame = lsq.theta_spread_deg; end
if TH_CONST && isfinite(lsq.theta_const)
  th_frame = th_geo_raw(:);
  q_frame = lsq.theta_const_q * double(isfinite(th_frame));
else
  [th_frame, q_frame] = H_smooth_profile(th_geo_raw, lsq.q_theta);
end
sub_frame = zw(:) + half > zv_frame;
th_frame(sub_frame) = NaN;
q_frame(sub_frame) = 0;

% --- segmentation in along-track distance
span = x_along(end) - x_along(1);
nseg = max(1, round(abs(span) / SEG_LEN));
fp = struct('zw', zw, 'th_frame', th_frame, 'q_frame', q_frame, ...
  'nseg', nseg, 'curved', curved, 'hspread', hspread, ...
  'ped_ant', ped_ant, 'track_az', track_az, 'lsq', lsq, ...
  'th_spread_frame', th_spread_frame);
if nseg == 1
  fp.th_seg = th_frame;
  fp.q_seg = q_frame;
  fp.dlam_seg = lsq.dlam;
  fp.resid_seg = lsq.resid;
  fp.seg_x = mean(x_along([1 end]));
  fp.seg_n = Nx;
  fp.th_spread_seg = th_spread_frame;   % the frame IS the one segment
  fp.held_seg = TH_CONST && isfinite(lsq.theta_const);
  fp.se_theta_seg = nan(Nw, 1); fp.se_dlam_seg = nan(Nw, 1);
  fp.jack_n = 0; fp.jack_edge = 0;
  fp.se_theta_n = zeros(Nw, 1); fp.se_theta_r = nan(Nw, 1);
  fp.se_theta_sat = false(Nw, 1);
  if JACK
    % The straight frame pass fitted the antenna-frame channels directly;
    % the replicates resample the same traces as geographic sub-block
    % moments (identical up to the track rotation on a straight line, and
    % exactly the curved path otherwise), around the frame's own axis.
    [~, Msub, nsub] = H_geo_moments(T, az_tr, 1:Nx, NBLK_ROT, CHAN);
  else
    Msub = {}; nsub = [];
  end
  % Same floor as the per-segment loop below: a delete-one jackknife needs
  % at least three sub-blocks to have a variance at all, and with fewer
  % ptt.quadpolJackknife either indexes an empty cell or spends a full
  % estimator pass on a degenerate replicate to return all-NaN. A short
  % frame reaches this branch now that the pipeline asks for a jackknife
  % unconditionally, so it must decline rather than fail.
  if JACK && numel(Msub) >= 3
    og = segbase;
    og.pedestal = H_ped_field(ped_ant, az_tr, 1:Nx, PSI_FIT);
    ref = lsq;
    ref.theta0 = th_geo_raw(:);
    if isfinite(lsq.theta_const), ref.theta_const = lsq.theta_const + deg2rad(track_az) * ~curved; end
    J = ptt.quadpolJackknife(Msub, nsub, z, og, ref);
    fp.se_theta_seg = J.se_theta; fp.se_dlam_seg = J.se_dlam;
    fp.jack_n = J.n; fp.jack_edge = J.n_edge;
    fp.se_theta_n = J.n_theta; fp.se_theta_r = J.r_theta;
    fp.se_theta_sat = J.sat_theta;
  end
  return
end

xb = linspace(min(x_along), max(x_along), nseg + 1);
th_raw = nan(Nw, nseg);
th_spread = nan(1, nseg);   % per-segment axis spread: is the assumption true here?
held_seg = false(1, nseg);  % segment fitted as ONE axis (constant mode, vote succeeded)
q_raw = zeros(Nw, nseg);
dlam_seg = nan(Nw, nseg);
resid_seg = nan(Nw, nseg);
seg_x = nan(1, nseg);
seg_n = zeros(1, nseg);
se_theta_seg = nan(Nw, nseg);
se_dlam_seg = nan(Nw, nseg);
se_theta_n = zeros(Nw, nseg);
se_theta_r = nan(Nw, nseg);
se_theta_sat = false(Nw, nseg);
jack_n = zeros(1, nseg);
jack_edge = zeros(1, nseg);
zv_seg = inf(1, nseg);
for s = 1:nseg
  js = find(x_along >= xb(s) & (x_along < xb(s+1) | s == nseg));
  seg_n(s) = numel(js);
  seg_x(s) = (xb(s) + xb(s+1)) / 2;
  zv_seg(s) = H_zvalid(Z_BED, js, BED_MARGIN, zv_frame);
  if numel(js) < 32, continue; end
  [Mg, Msub, nsub] = H_geo_moments(T, az_tr, js, NBLK_ROT, CHAN);
  if isempty(Mg), continue; end
  os = segbase;   % carries theta_const: one axis per SEGMENT, constant in depth
  os.z_valid = zv_seg(s);
  os.pedestal = H_ped_field(ped_ant, az_tr, js, PSI_FIT);
  o = ptt.quadpolFabricLS(struct('M', Mg), z, os);
  if nnz(isfinite(o.theta0)) < MIN_SEG_W, continue; end
  th_raw(:, s) = o.theta0;
  if isfield(o, 'theta_spread_deg'), th_spread(s) = o.theta_spread_deg; end
  held_seg(s) = TH_CONST && isfinite(o.theta_const);
  if held_seg(s)
    % one axis, one weight (the pooled contrast) - see the frame profile
    q_raw(:, s) = o.theta_const_q * double(isfinite(o.theta0));
  else
    qs = o.q_theta;
    qs(~isfinite(qs) | qs < 0) = 0;
    qs(~isfinite(o.theta0)) = 0;
    q_raw(:, s) = qs;
  end
  dlam_seg(:, s) = o.dlam;
  resid_seg(:, s) = o.resid;
  if JACK && numel(Msub) >= 3
    J = ptt.quadpolJackknife(Msub, nsub, z, os, o);
    se_theta_seg(:, s) = J.se_theta;
    se_dlam_seg(:, s) = J.se_dlam;
    se_theta_n(:, s) = J.n_theta;
    se_theta_r(:, s) = J.r_theta;
    se_theta_sat(:, s) = J.sat_theta;
    jack_n(s) = J.n;
    jack_edge(s) = J.n_edge;
  end
end

% --- smooth each segment in depth (as the frame profile is) on the
% q-weighted doubled-angle phasor; never unwrap. There is deliberately NO
% smoothing ACROSS segments: a ~2 km segment already pools thousands of
% traces, and any lateral kernel drags a pure segment toward its
% neighbour - measured 7 deg of pull on a 45 deg margin step, exactly
% where lateral resolution is the point. Continuity in the handoff comes
% from the weighted phasor interpolation in ptt.thetaProfileAt.
th_seg = nan(Nw, nseg);
q_seg = zeros(Nw, nseg);
for s = 1:nseg
  if held_seg(s)
    th_seg(:, s) = th_raw(:, s);   % constant by construction; nothing to smooth
    q_seg(:, s) = q_raw(:, s);
  else
    [th_seg(:, s), q_seg(:, s)] = H_smooth_profile(th_raw(:, s), q_raw(:, s));
  end
end
% windows at or below the segment's own z_valid are not ice: no axis, no
% weight, and excluded from both fallbacks below
sub_seg = zw(:) + half > zv_seg;
th_seg(sub_seg) = NaN;
q_seg(sub_seg) = 0;
% fallback, two levels, value-only (weight stays 0 so the handoff
% interpolation cannot be dragged): a dead segment window ABOVE its
% z_valid first BORROWS the q-weighted phasor mean of the segments that
% are alive at that window - the axis field is the slowly-varying
% quantity, and a live neighbour beats any dead pooled fit - and only
% where no segment is alive does the frame profile fill in
dead = ~(q_seg > 0);
ph_mat = q_seg .* exp(2i*th_seg);
ph_mat(dead | ~isfinite(ph_mat)) = 0;
ph_row = sum(ph_mat, 2);
w_row = sum(q_seg .* ~dead, 2);
borrow = dead & ~sub_seg & repmat(w_row > 0, 1, nseg);
thb = repmat(0.5 * angle(ph_row), 1, nseg);
th_seg(borrow) = thb(borrow);
fb = dead & ~borrow & ~sub_seg & repmat(isfinite(th_frame(:)), 1, nseg);
thf = repmat(th_frame(:), 1, nseg);
th_seg(fb) = thf(fb);

fp.th_seg = th_seg;
fp.q_seg = q_seg;
fp.dlam_seg = dlam_seg;
fp.resid_seg = resid_seg;
fp.seg_x = seg_x;
fp.seg_n = seg_n;
fp.th_spread_seg = th_spread;
fp.held_seg = held_seg;
fp.se_theta_seg = se_theta_seg;
fp.se_dlam_seg = se_dlam_seg;
fp.se_theta_n = se_theta_n;
fp.se_theta_r = se_theta_r;
fp.se_theta_sat = se_theta_sat;
fp.jack_n = jack_n;
fp.jack_edge = jack_edge;
end

% -------------------------------------------------------------------------
function v = H_opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end

function [ths, qs] = H_smooth_profile(th_raw, q)
% q-weighted doubled-angle phasor smoothing over 5 windows in depth; the
% same kernel and threshold as the unsegmented pipeline handoff.
ok = isfinite(th_raw);
qw = q;
qw(~isfinite(qw) | qw < 0) = 0;
qw(~ok) = 0;
ph = qw .* exp(2i*th_raw);
ph(~ok) = 0;
ks = ones(5, 1);
num_p = conv(ph, ks, 'same');
den_p = conv(qw, ks, 'same');
ths = nan(size(th_raw));
qs = zeros(size(th_raw));
oks = den_p > 0.25 & abs(num_p) > 0;
ths(oks) = 0.5 * angle(num_p(oks));
qs(oks) = den_p(oks);
if nnz(oks) < 2 && nnz(ok) >= 2
  % too weak to smooth but not empty: pass the raw profile through, as
  % the unsegmented handoff does
  ths(ok) = th_raw(ok);
  qs(ok) = max(qw(ok), eps);
end
end

function [Mg, Msub, nsub] = H_geo_moments(T, az_tr, js, nblk_rot, chan)
% pooled geographic-frame moments over the trace subset js: heading-
% rotated per sub-block, trace-count weighted; sub-blocks under 16 traces
% are skipped (too few for a moment average worth rotating). The
% sub-block moments and their trace counts are returned too - they are
% the jackknife's resampling units, and Mg is exactly their weighted mean.
% The counts are PER DEPTH, [Nt x n]: where traces are blanked below the
% bed (ptt.maskBelowBed) a sub-block contributes only the traces still in
% ice at each depth, and a depth with none carries no weight rather than
% turning the pooled moment NaN. Unmasked, every column is constant and
% the pooling is the plain trace-count mean it always was.
Mg = [];
wsum = [];
Msub = {};
nsub = [];
for jr = 1:nblk_rot:numel(js)
  jr1 = min(jr + nblk_rot - 1, numel(js));
  if jr1 - jr < 16, continue; end
  jj = js(jr:jr1);
  Sb = struct();
  for k = 1:4, Sb.(chan{k}) = T.(chan{k})(:, jj); end
  Mb = ptt.quadpolMoments(Sb, [1 numel(jj)]);
  hb = mod(rad2deg(angle(mean(exp(2i*deg2rad(az_tr(jj))))))/2, 180);
  Mr = ptt.rotateMoments(Mb, -deg2rad(hb));
  cnt = sum(isfinite(Sb.(chan{1})), 2);        % traces in ice, per depth
  Mz = Mr;
  Mz(~isfinite(Mz)) = 0;
  if isempty(Mg), Mg = zeros(size(Mr)); wsum = zeros(size(cnt)); end
  Mg = Mg + Mz .* cnt;
  wsum = wsum + cnt;
  Msub{end+1} = Mr; %#ok<AGROW>
  nsub(:, end+1) = cnt; %#ok<AGROW>
end
if ~isempty(wsum)
  Mg = Mg ./ wsum;                             % 0/0 -> NaN where no ice
end
end

function pfld = H_ped_field(ped_ant, az_tr, js, psi_fit)
% the antenna-frame pedestal mixed over the heading distribution of the
% trace subset: rotation carries sin2(phi) to sin2(psi - h)
mh2 = mean(exp(2i*deg2rad(az_tr(js))));
mh4 = mean(exp(4i*deg2rad(az_tr(js))));
c2h = real(mh2); s2h = imag(mh2);
c4h = real(mh4); s4h = imag(mh4);
pfld = (ped_ant(1) + 1i*ped_ant(2)) ...
  * (c2h*sin(2*psi_fit(:)) - s2h*cos(2*psi_fit(:))) ...
  + ped_ant(3) * (0.5 - 0.5*(c4h*cos(4*psi_fit(:)) + s4h*sin(4*psi_fit(:))));
end

function zv = H_zvalid(z_bed, js, margin, fallback)
% The deepest window bottom at which at least half of traces js are still in
% ice: their MEDIAN picked bed, less the margin. Unpicked traces count as
% ice (they are not masked either); with no pick at all the caller's
% fallback stands (Inf when no bed is known).
zv = fallback;
if isempty(z_bed), return; end
b = z_bed(js);
if ~any(isfinite(b)), return; end
b(~isfinite(b)) = Inf;
zv = median(b) - margin;
end
