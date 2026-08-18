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
% weighted phasor per segment; weightless (fallback) segments enter with a
% small epsilon so a fully-dead row still interpolates its frame values
w = fp.q_seg;
w(~isfinite(fp.th_seg)) = 0;
w = w + 1e-6 * double(isfinite(fp.th_seg));
ph = w .* exp(2i*fp.th_seg);
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
