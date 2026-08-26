function tf = pedestalFailed(ped)
%PEDESTALFAILED True when a pedestal estimate is not a measurement.
%
% tf = ptt.pedestalFailed(ped)
%
% The one predicate every consumer of a frame pedestal tests, so that the
% writer's marker and the reader's test cannot drift apart. It is true for
% a non-finite estimate AND for ptt.pedestalMarker() - the exact zeros a
% failed frame-pass fit is recorded as - because both mean the same thing:
% this frame has no pedestal, so anything that would otherwise be anchored
% to it must fall back instead (in the pipeline, to each block fitting its
% own). Anything not shaped like a [1 x 3] coefficient triple is failed
% too, since it cannot be used as one.
%
% See also ptt.pedestalMarker, ptt.quadpolFrameTheta, ptt.quadpolFabricLS.

p = ped(:).';
tf = numel(p) ~= 3 || ~all(isfinite(p)) || isequal(p, ptt.pedestalMarker());

end
