function prof = thetaProfileAt(fp, x)
%THETAPROFILEAT Geographic theta0 handoff profile at an along-track position.
%
% prof = ptt.thetaProfileAt(fp, x)
%
% Interpolates the segmented frame pass (ptt.quadpolFrameTheta) at
% along-track position x [m], returning the struct('z', 'theta') form
% that ptt.quadpolFabricLS accepts as opts.theta0 (the caller subtracts
% the block heading to reach the block's antenna frame). Returns [] when
% the frame pass produced nothing usable, which tells the block fit to
% estimate its own axis - the same convention the unsegmented pipeline
% used.
%
% Interpolation is across segment centres on the DOUBLED-ANGLE PHASOR,
% weighted by the segment weights so a fallback-filled segment (weight 0)
% between live neighbours is bridged by them rather than pulling the
% profile toward the frame average; where every segment is weightless the
% plain phasor of the fallback values (the frame profile) survives.
% Positions beyond the end segments clamp, matching the depth-clamping
% rule in the estimator: a linear extrapolation of a phasor can pass
% near zero and land on an arbitrary angle.
%
% BELOW A SEGMENT'S BED. The frame pass leaves a segment's windows at or
% below its own median bed (its z_valid) empty - no axis, no weight - so
% the saved profiles never show fabric where there is no ice. The handoff
% treats the two modes differently there. A HELD segment asserts ONE axis
% for its column, so it is handed off at every depth with its one weight:
% a block between two held segments interpolates two constants and never
% switches axis at the shallower neighbour's bed (it used to, taking the
% deeper neighbour's axis alone between the two beds). A FREE segment
% says nothing about its axis below its bed, so there a block's rows are
% bridged from the live neighbours' deeper ice, the convention for any
% dead row; the block's own per-trace bed mask (ptt.maskBelowBed), not
% this profile, is what keeps the block's sub-bed windows out.

if isempty(fp) || ~isfield(fp, 'th_seg') || isempty(fp.th_seg)
  prof = [];
  return
end
if fp.nseg == 1
  ok = isfinite(fp.th_seg(:, 1));
  if nnz(ok) < 2, prof = []; return; end
  prof = struct('z', fp.zw(ok), 'theta', fp.th_seg(ok, 1));
  return
end

xq = min(max(x, fp.seg_x(1)), fp.seg_x(end));
% a held segment's one axis and weight at every depth (see BELOW A
% SEGMENT'S BED)
th = fp.th_seg;
w = fp.q_seg;
if isfield(fp, 'held_seg')
  for s = find(fp.held_seg(:).')
    live = isfinite(th(:, s)) & w(:, s) > 0;
    if ~any(live), continue; end
    % the live windows already carry the axis; only the empty ones change
    th(~live, s) = 0.5 * angle(sum(w(live, s) .* exp(2i*th(live, s))));
    w(~live, s) = max(w(live, s));
  end
end
% weighted phasor per segment; weightless (fallback) segments enter with a
% small epsilon so a fully-dead row still interpolates its frame values
w(~isfinite(th)) = 0;
w = w + 1e-6 * double(isfinite(th));
ph = w .* exp(2i*th);
ph(~isfinite(ph)) = 0;
num = interp1(fp.seg_x(:), ph.', xq, 'linear').';
den = interp1(fp.seg_x(:), w.', xq, 'linear').';
ok = den > 0 & abs(num) > 0 & isfinite(num);
if nnz(ok) < 2, prof = []; return; end
zp = fp.zw(ok);
tp = 0.5 * angle(num(ok));
wp = den(ok);
% Robust end rows: the estimator clamps every window beyond the profile's
% range to the END value, so when the profile dies early (deep coherence
% loss) a single noisy terminal row would be propagated to every deeper
% window of every block. Replace each end value with the q-weighted
% phasor mean of the outermost min(5, n) rows.
K = min(5, numel(tp));
t_bot = 0.5 * angle(sum(wp(end-K+1:end) .* exp(2i*tp(end-K+1:end))));
t_top = 0.5 * angle(sum(wp(1:K) .* exp(2i*tp(1:K))));
tp(end) = t_bot;
tp(1) = t_top;
prof = struct('z', zp, 'theta', tp);
end
