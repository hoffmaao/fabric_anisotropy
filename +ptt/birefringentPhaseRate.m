function gpd1 = birefringentPhaseRate(fc, eps_perp, deps)
%BIREFRINGENTPHASERATE One-way relative phase per metre per unit dlam.
%
% gpd1 = ptt.birefringentPhaseRate(fc, eps_perp, deps)
%
% The x-vs-y phase a birefringent layer accumulates in ONE pass, per metre
% and per unit horizontal eigenvalue difference dlam:
%
%   gpd1 = pi * fc * Deps' / (sqrt(eps') * c)
%
% from Dn = Deps' dlam / (2 sqrt(eps')) and a one-way phase 2 pi fc Dn / c.
% A down-and-up sandwich doubles it, which is the two-way gradient
% ptt.ershadiFabric inverts in its eq. (9)/(10) step - the two are exact
% reciprocals of one another, so they must be computed from the same
% constants or a retrieval and the model it is compared against land on
% different depths.
%
% SINGLE OWNER. ptt.fujitaModel propagates phase with this, and
% ptt.ershadiInverse locates the anti-phase depths for its eq. (13)
% cross-check with it. Those two must agree exactly or the cross-check is
% evaluated at depths the model never puts a node at, which is the same
% class of silent desync that ptt.ershadiInverse's H_const guards for
% fc/eps_perp/deps themselves - so the constant they feed lives here rather
% than being written out at each call site.
%
% Inputs
%   fc        centre frequency (Hz)
%   eps_perp  perpendicular permittivity (the paper's 3.15, not eps_bar)
%   deps      dielectric anisotropy (ptt.constants deps, 0.034)
%
% See also ptt.fujitaModel, ptt.ershadiInverse, ptt.ershadiFabric.

% Three positional scalars of the same type invite a transposed call, and a
% transposed call here is silent - it returns a plausible number and shifts
% every modelled phase. eps_perp is a permittivity (> 1) and deps a small
% anisotropy (< 1), so the common swap cannot survive these checks.
if ~isscalar(fc) || ~isfinite(fc) || fc <= 0
  error('ptt:birefringentPhaseRate:fc', 'fc must be a positive scalar Hz');
end
if ~isscalar(eps_perp) || ~isfinite(eps_perp) || eps_perp <= 1
  error('ptt:birefringentPhaseRate:epsPerp', ...
    'eps_perp must be a scalar permittivity > 1, got %g', eps_perp);
end
if ~isscalar(deps) || ~isfinite(deps) || deps <= 0 || deps >= 1
  error('ptt:birefringentPhaseRate:deps', ...
    'deps must be a scalar anisotropy in (0, 1), got %g', deps);
end

C0 = ptt.constants();
gpd1 = pi * fc * deps / (sqrt(eps_perp) * C0.c * 1e9);

end
