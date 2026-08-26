function p = pedestalMarker()
%PEDESTALMARKER The value saved in place of a pedestal that did not converge.
%
% p = ptt.pedestalMarker()
%
% Exactly [0 0 0]. A converged fit never lands on all three coefficients
% being exactly zero, so this value is unambiguous as an audit marker: a
% frame whose saved ls_pedestal reads it is a frame with NO pedestal
% measurement, and every consumer must drop it rather than use it (see
% docs/method.md). The value is fixed - products and audits already key on
% it - so it lives here once instead of as a literal in each writer, and
% ptt.pedestalFailed is the matching test.
%
% See also ptt.pedestalFailed, ptt.quadpolFrameTheta.

p = [0 0 0];

end
