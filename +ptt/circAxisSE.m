function [s, m, Rbar, sat] = circAxisSE(th, min_rep)
%CIRCAXISSE Delete-one jackknife standard error of an AXIS, bounded.
%
% [s, m, Rbar, sat] = ptt.circAxisSE(th, min_rep)
%
% th       replicate axis angles [rad], any shape; NaNs ignored. Each is an
%          AXIS, unique modulo pi, so all arithmetic here is circular on
%          the DOUBLED angle and nothing is ever unwrapped.
% min_rep  fewer than this many finite replicates returns NaN (default 3)
%
% s        standard error of the axis [rad], or NaN when the replicates
%          carry no axis the statistic can resolve (see ABSTENTION)
% m        replicates that contributed
% Rbar     resultant length of the doubled angles, in [0, 1]; NaN when
%          there were too few replicates to form one
% sat      true when s is NaN BECAUSE the replicates disagreed past what
%          this statistic can resolve, false otherwise. With m and Rbar it
%          separates the three things a NaN can mean: too few replicates
%          (m < min_rep, sat false), replicates that scattered beyond the
%          resolvable range (sat true), and - when s is finite - a genuine
%          measurement. They are different statements about the ice and
%          the instrument and must not collapse into one bare NaN.
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
% still assert a measurement that was not made, and a bounded statistic
% stops discriminating well before its bound. MEASURED over 3000 draws per
% point: at m = 22 the reported value sits at 36.7 deg for every true
% scatter from 16 deg to 40 deg, and at m = 11 it plateaus near 35 deg
% from 16 deg upward - a factor of two of real uncertainty compressed into
% a number that still reads as a measurement. So there are two abstentions,
% and both set sat:
%
%  - the axes themselves are indistinguishable from a uniform spread, by
%    Rayleigh on the doubled angle (m*Rbar^2 is ~1 for uniform samples)
%    against its 5% point, 2.996;
%  - the STRETCHED resultant has decayed past the point where it still
%    tracks the scatter, max(0.4058, 1/sqrt(m)). Two effects set that, and
%    the binding one wins. RESOLUTION: with wrapped-normal deviations
%    Rj = exp(-x^2/2) for x = sqrt(m-1)*sigma, so ds/dsigma falls away
%    from its small-scatter value as x grows; requiring at least half that
%    sensitivity gives x = 1.343, Rj = 0.4058, s = 31.2 deg, and beyond it
%    a reported number no longer distinguishes 16 deg of scatter from 40.
%    SAMPLING: m phasors in random directions give E[Rj^2] = 1/m, so at
%    small m an Rj below 1/sqrt(m) is what noise alone produces. The floor
%    is therefore m-dependent, as the plateau is - sampling binds up to
%    m = 6, resolution from m = 7 on - and caps what can be reported at
%    30.1 deg for m = 5 and 31.2 deg above it.
%
% Widening the bound to avoid saturating would be the wrong fix. Saturation
% is correct behaviour for a bounded statistic on a circular quantity; the
% defect would be reporting a value from inside the band where the bound
% has already destroyed the information.
%
% See also ptt.quadpolJackknife.

if nargin < 2 || isempty(min_rep), min_rep = 3; end

RAYLEIGH_05 = 2.9957;   % -log(0.05), the 5% point of m*Rbar^2
RESOLVE_MIN = 0.4058;   % stretched resultant at half sensitivity

ok = isfinite(th);
m = nnz(ok);
s = NaN;
Rbar = NaN;
sat = false;
if m < min_rep, return; end

u = 2 * th(ok);
u = u(:);
zb = mean(exp(1i * u));
Rbar = abs(zb);

if m * Rbar^2 < RAYLEIGH_05, sat = true; return; end

d = angle(exp(1i * (u - angle(zb))));        % doubled-angle deviations
Rj = abs(mean(exp(1i * sqrt(m - 1) * d)));
if Rj <= max(RESOLVE_MIN, 1 / sqrt(m)), sat = true; return; end
s = sqrt(2 * (1 - Rj)) / 2;
end
