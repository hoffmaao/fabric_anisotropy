function out = strainRateAxes(in, opts)
%STRAINRATEAXES Surface principal strain-rate axes from a velocity field.
%
% out = ptt.strainRateAxes(in, opts)
%
% Returns the principal axes and rates of the horizontal surface
% strain-rate tensor, for referencing fabric orientation against the
% quantity that actually orients fabric.
%
% WHY STRAIN RATE AND NOT FLOW DIRECTION. Fabric is built by the
% deformation history, not by the instantaneous velocity vector. The two
% coincide where flow is fast and straight, and part company exactly where
% this project works hardest: at a divide the speed is small and its
% direction poorly determined, while the strain-rate axes stay well
% defined; in a shear margin flow is nearly uniaxial while the strain is
% dominated by simple shear whose principal axes sit at 45 degrees to it.
% A flow-referenced angle in a margin is therefore off by 45 degrees from
% the deformation the ice actually felt.
%
% in is a struct, either
%   .exx .eyy .exy   strain-rate components [1/a or 1/s], any consistent
%                    unit; only ratios and axes are used here
% or
%   .vx .vy          velocity components on the projection GRID [m/a],
%                    2-D arrays
%   .dx .dy          grid spacing [m] (scalars; dy defaults to dx)
%
% Gradients are taken with respect to grid x (east) and grid y (north) of
% the projection, so every azimuth returned is measured from GRID north,
% clockwise - the same convention as atan2(vx, vy). Grid north is not true
% north at these latitudes; ptt.orientationToStrain refuses to mix the two
% frames rather than silently differencing them.
%
% opts:
%   .smooth_px (0)  boxcar half-width in pixels applied to vx, vy before
%                   differencing. Differencing amplifies noise, and a
%                   velocity mosaic is noisy; 0 leaves the field alone.
%
% out.az_ext    azimuth of the MAXIMUM principal (most extensional) axis
%               [rad from grid north], modulo pi - it is an axis
% out.az_cmp    azimuth of the minimum principal axis, always az_ext+pi/2
% out.e_ext     maximum principal strain rate (algebraically largest)
% out.e_cmp     minimum principal strain rate
% out.e_eff     effective strain rate, sqrt(0.5*(exx^2+eyy^2)+exy^2)
% out.regime    'extension' | 'compression' | 'shear', from the signs and
%               ratio of the principal rates
% out.shear_frac  |e_ext+e_cmp| / (|e_ext|+|e_cmp|) subtracted from 1:
%               0 for pure uniaxial, 1 for pure shear (equal and opposite)
%
% See also ptt.orientationToStrain, ptt.orientationToFlow.
if nargin < 2, opts = struct(); end
if isfield(in, 'exx') && isfield(in, 'eyy') && isfield(in, 'exy')
  exx = double(in.exx); eyy = double(in.eyy); exy = double(in.exy);
elseif isfield(in, 'vx') && isfield(in, 'vy')
  if ~isfield(in, 'dx') || isempty(in.dx)
    error('ptt:strainRateAxes:spacing', 'supply in.dx (grid spacing, m)');
  end
  dx = double(in.dx);
  if isfield(in, 'dy') && ~isempty(in.dy), dy = double(in.dy); else, dy = dx; end
  vx = double(in.vx); vy = double(in.vy);
  sm = H_opt(opts, 'smooth_px', 0);
  if sm > 0, vx = H_box(vx, sm); vy = H_box(vy, sm); end
  % gradient() returns d/dcol then d/drow; columns are grid x (east) and
  % rows are grid y (north) for a north-up raster, so the row derivative
  % is taken against +y explicitly rather than by whichever way the
  % raster happens to be stored.
  [dvx_dx, dvx_dy] = gradient(vx, dx, dy);
  [dvy_dx, dvy_dy] = gradient(vy, dx, dy);
  exx = dvx_dx; eyy = dvy_dy; exy = 0.5 * (dvx_dy + dvy_dx);
else
  error('ptt:strainRateAxes:input', ...
    'in needs either .exx/.eyy/.exy or .vx/.vy with .dx');
end

mn = (exx + eyy) / 2;
df = (exx - eyy) / 2;
rr = hypot(df, exy);
e1 = mn + rr;            % algebraically largest principal rate
e2 = mn - rr;

% Principal direction measured from +x (grid east), counterclockwise:
%   phi = 0.5 * atan2(2*exy, exx - eyy)
% converted to an azimuth from grid north, clockwise, and reduced modulo
% pi because a principal direction is an AXIS, not a direction.
phi = 0.5 * atan2(2*exy, exx - eyy);
az_ext = mod(pi/2 - phi, pi);

out = struct();
out.az_ext = az_ext;
out.az_cmp = mod(az_ext + pi/2, pi);
out.e_ext = e1;
out.e_cmp = e2;
out.e_eff = sqrt(0.5*(exx.^2 + eyy.^2) + exy.^2);
out.exx = exx; out.eyy = eyy; out.exy = exy;
den = abs(e1) + abs(e2);
out.shear_frac = 1 - abs(e1 + e2) ./ max(den, eps);
reg = repmat({'shear'}, size(e1));
reg(e1 > 0 & e2 >= 0) = {'extension'};
reg(e1 <= 0 & e2 < 0) = {'compression'};
if isscalar(reg), out.regime = reg{1}; else, out.regime = reg; end
end

function y = H_box(x, h)
k = ones(2*h+1, 1) / (2*h+1);
y = conv2(k, k.', x, 'same');
e = conv2(k, k.', ones(size(x)), 'same');
y = y ./ max(e, eps);
end

function v = H_opt(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
