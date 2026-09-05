function [s, m, Rbar] = circAxisSE(th, min_rep)
%CIRCAXISSE Delete-one jackknife standard error of an AXIS, bounded.
%
% [s, m, Rbar] = ptt.circAxisSE(th, min_rep)
%
% th       replicate axis angles [rad], any shape; NaNs ignored. Each is an
%          AXIS, unique modulo pi, so all arithmetic here is circular on
%          the DOUBLED angle and nothing is ever unwrapped.
% min_rep  fewer than this many finite replicates returns NaN (default 3)
%
% s        standard error of the axis [rad], or NaN when the replicates
%          carry no axis at all (see ABSTENTION)
% m        replicates that contributed, so a NaN from "too few" can be told
%          apart from a NaN from "they did not agree"
% Rbar     resultant length of the doubled angles, in [0, 1]
%
% WHY NOT THE LINEAR JACKKNIFE. The delete-one formula
% sqrt((m-1)/m * sum(d.^2)) is right for an estimate living on an unbounded
% line. An axis does not: |d| <= pi/2 by construction, so that sum has a
% ceiling of (pi/2)*sqrt(m-1) - 180 deg at m = 5, 285 deg at m = 11, 412
% deg at m = 22. It reports MORE than a half turn of uncertainty for a
% quantity defined modulo a half turn, and it does so at replicate scatter
% barely above (m = 11) or even below (m = 22) the circular-uniform level.
% MEASURED on the 35 constant-orientation products on disk: median 79 deg,
% 64% of segments above the 52 deg uniform SD, 24% above 180 deg, maximum
% 1290 deg, while the same segments' full-data fits were healthy.
%
% WHAT THIS COMPUTES INSTEAD. The dispersion comes from a RESULTANT, which
% is bounded in [0, 1] by construction, so the statistic cannot exceed its
% bound in the first place rather than being clamped at one. Each deviation
% from the mean direction is stretched by the delete-one factor sqrt(m-1)
% and the resultant of the STRETCHED deviations is taken:
%
%   Rj = |mean(exp(1i * sqrt(m-1) * d))|,   s = sqrt(2*(1 - Rj)) / 2
%
% The trailing /2 carries the doubled angle back to the axis. For tightly
% clustered replicates 1 - Rj -> mean(d.^2)/2 and s -> sqrt(m-1)*rms(d)/2,
% which is the linear formula exactly, so the delete-one inflation is kept
% and nothing is made optimistic. As the replicates spread, the stretched
% phasors wrap and Rj falls to 0, so s saturates at sqrt(2)/2 = 40.5 deg -
% under the 52 deg uniform SD at every m, with no sqrt(m-1) growth left.
%
% ABSTENTION. Saturating is not enough: a value near the ceiling would
% still assert a measurement that was not made. Replicates that are
% indistinguishable from a uniform spread of axes determine no axis, so
% they return NaN, following the convention used throughout this package.
% The test is Rayleigh's on the doubled angle - m*Rbar^2 is ~1 for uniform
% samples - against its 5% point, 2.996.
%
% See also ptt.quadpolJackknife.

if nargin < 2 || isempty(min_rep), min_rep = 3; end

RAYLEIGH_05 = 2.9957;   % -log(0.05), the 5% point of m*Rbar^2

ok = isfinite(th);
m = nnz(ok);
s = NaN;
Rbar = NaN;
if m < min_rep, return; end

u = 2 * th(ok);
u = u(:);
zb = mean(exp(1i * u));
Rbar = abs(zb);

if m * Rbar^2 < RAYLEIGH_05, return; end

d = angle(exp(1i * (u - angle(zb))));        % doubled-angle deviations
Rj = abs(mean(exp(1i * sqrt(m - 1) * d)));
s = sqrt(2 * max(1 - Rj, 0)) / 2;
end
