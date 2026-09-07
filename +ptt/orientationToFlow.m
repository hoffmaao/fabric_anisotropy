function out = orientationToFlow(theta0, flow, opts)
%ORIENTATIONTOFLOW Fabric axis expressed relative to the local flow direction.
%
% out = ptt.orientationToFlow(theta0, flow, opts)
%
% Converts a fabric AXIS into the angle it makes with the local surface
% flow, which is how the NEGIS orientations are reported (against MEaSUREs
% multi-year velocity, with GNSS on the antenna fixing the radar's own
% heading). The reported trend there is toward flow alignment with
% increasing upstream distance from EastGRIP, so the sign and range
% conventions below decide whether that trend is visible at all.
%
% AN AXIS IS NOT A DIRECTION, and this is the whole subtlety. theta0 is
% defined modulo 180 degrees: it has no head or tail, so "aligned with
% flow" and "anti-aligned" are the same statement. Flow is a true
% direction modulo 360. The angle between them is therefore confined to
% [0, 90] degrees - 0 is flow-parallel, 90 is flow-perpendicular - and any
% value outside that range means a convention was mixed up somewhere
% upstream. Everything here is computed on the DOUBLED angle and never by
% unwrapping.
%
% GRID NORTH IS NOT TRUE NORTH, and at these latitudes the difference is
% large. MEaSUREs velocities are components on a polar-stereographic grid
% (EPSG:3413 in Greenland, 3031 in Antarctica), so an azimuth built from
% atan2(vx, vy) is measured from GRID north. A theta0 in true-north
% azimuth cannot be differenced against it without the convergence
% correction. Rather than guess, this function requires the caller to say
% which frame each input is in and refuses to mix them silently.
%
% theta0  fabric axis [rad], modulo pi, scalar or array
% flow    either a scalar/array flow azimuth [rad] in the SAME frame, or a
%         struct with fields .vx and .vy (velocity components on the
%         projection grid), from which a grid azimuth is formed
% opts:
%   .theta_frame  'true' | 'grid'   frame of theta0   (REQUIRED)
%   .flow_frame   'true' | 'grid'   frame of flow     (REQUIRED when flow
%                                   is an azimuth; forced 'grid' for vx/vy)
%   .convergence  [rad] grid-to-true convergence at the site, to be added
%                 to a grid azimuth to obtain a true one. Required only
%                 when the two frames differ.
%   .v_min        (0) speed below which the flow direction is treated as
%                 undefined and the result abstains (NaN). A direction
%                 from a near-zero velocity is noise, not slow flow.
%
% out.angle_deg   angle between axis and flow, in [0, 90], NaN where
%                 abstained
% out.flow_az     the flow azimuth actually used [rad]
% out.aligned     logical, angle < 45 (closer to along-flow than across)
% out.speed       flow speed where vx/vy were supplied, else []
%
% See also ptt.nymandTwoStep.
if nargin < 3, opts = struct(); end
tf = H_req(opts, 'theta_frame', 'ptt:orientationToFlow:frame', ...
  'opts.theta_frame must be ''true'' or ''grid'' - see the header on grid vs true north');

speed = [];
if isstruct(flow)
  if ~isfield(flow, 'vx') || ~isfield(flow, 'vy')
    error('ptt:orientationToFlow:vxvy', 'flow struct needs .vx and .vy');
  end
  vx = double(flow.vx); vy = double(flow.vy);
  speed = hypot(vx, vy);
  flow_az = atan2(vx, vy);          % azimuth from GRID north, clockwise
  ff = 'grid';
else
  flow_az = double(flow);
  ff = H_req(opts, 'flow_frame', 'ptt:orientationToFlow:frame', ...
    'opts.flow_frame must be ''true'' or ''grid'' when flow is an azimuth');
end

if ~strcmp(tf, ff)
  if ~isfield(opts, 'convergence') || isempty(opts.convergence)
    error('ptt:orientationToFlow:convergence', ...
      ['theta0 is in the %s frame and flow is in the %s frame; supply ' ...
       'opts.convergence (grid-to-true, rad) rather than differencing them'], tf, ff);
  end
  cvg = double(opts.convergence);
  if strcmp(ff, 'grid'), flow_az = flow_az + cvg; else, flow_az = flow_az - cvg; end
end

v_min = H_opt(opts, 'v_min', 0);
th = double(theta0);
% Angle between an AXIS and a DIRECTION, on the doubled angle so the
% axis's 180-degree ambiguity is handled by construction. The doubled
% difference folds head and tail onto the same point; halving it and
% taking the magnitude lands in [0, 90] with no unwrapping anywhere.
d = angle(exp(2i * (th - flow_az))) / 2;
ang = abs(rad2deg(d));

if ~isempty(speed) && v_min > 0
  ang(speed < v_min) = NaN;
end

out = struct();
out.angle_deg = ang;
out.flow_az = flow_az;
out.aligned = ang < 45;
out.speed = speed;
end

function v = H_req(s, f, id, msg)
if ~isstruct(s) || ~isfield(s, f) || isempty(s.(f)), error(id, '%s', msg); end
v = lower(char(s.(f)));
if ~ismember(v, {'true', 'grid'}), error(id, '%s', msg); end
end

function v = H_opt(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
