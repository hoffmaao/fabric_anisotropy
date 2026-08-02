function [depth, twtt] = twttDepthMap(par, npts)
%TWTTDEPTHMAP Vertical two-way traveltime versus depth below the surface.
%   [depth, twtt] = TWTTDEPTHMAP(par, npts) integrates the y-eigenslowness
%   of the firn-ice column model par (see ptt.defaultParams) from the
%   surface downward on npts points (default 4001). depth is in meters
%   below the surface, twtt in SECONDS (two-way). Birefringence is far too
%   weak to matter for this mapping.

if nargin < 2 || isempty(npts)
  npts = 4001;
end

zhat  = linspace(1, 0, npts).';
P     = ptt.columnProfiles(par, zhat);
depth = par.H * (1 - zhat);
twtt  = 2e-9 * cumtrapz(depth, P.S(:, 2)); % ptt slowness is ns/m

end
