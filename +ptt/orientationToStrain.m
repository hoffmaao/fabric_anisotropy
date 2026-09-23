function out = orientationToStrain(theta0, strain, opts)
%ORIENTATIONTOSTRAIN Fabric axis referenced to the principal strain-rate axes.
%
% out = ptt.orientationToStrain(theta0, strain, opts)
%
% Reports the fabric axis against the surface principal strain-rate axes
% rather than against the flow direction, because strain is what orients
% fabric. Where flow is fast and straight the two references agree; at a
% divide the flow direction is poorly determined while the strain axes are
% not, and in a shear margin the principal axes sit 45 degrees off the
% flow, so a flow-referenced angle there is wrong by that amount.
%
% BOTH INPUTS ARE AXES. A fabric axis is modulo 180 degrees and so is a
% principal strain direction, so the angle between them lies in [0, 90]:
% 0 means the fabric axis is parallel to the reference principal axis and
% 90 means perpendicular. The fold is done on the doubled angle by the
% shared private helper, never by unwrapping.
%
% WHICH PRINCIPAL AXIS IS THE REFERENCE, and why it is stated rather than
% assumed: opts.reference selects 'extension' (default, the algebraically
% largest principal rate) or 'compression'. The two answers are exact
% complements, angle_cmp = 90 - angle_ext, so reporting one without saying
% which is an invitation to read a perpendicular result as a parallel one.
% The output carries the choice back in out.reference for that reason.
%
% NO PHYSICAL SIGN IS IMPOSED. This function does not assert whether ice
% fabric should end up parallel or perpendicular to extension; that
% depends on the deformation regime and on which eigenvalue of the
% orientation tensor dlam names. It reports the geometry and leaves the
% interpretation to the caller.
%
% theta0  fabric axis [rad], modulo pi, scalar or array
% strain  either a struct from ptt.strainRateAxes, or an azimuth [rad] of
%         the reference principal axis directly
% opts:
%   .theta_frame  'true' | 'grid'  frame of theta0 (REQUIRED)
%   .strain_frame 'true' | 'grid'  frame of the strain axes; defaults to
%                 'grid', because ptt.strainRateAxes differentiates on the
%                 projection grid and returns grid azimuths
%   .convergence  [rad] grid-to-true convergence, required only when the
%                 two frames differ
%   .reference    'extension' (default) | 'compression'
%   .e_min        (0) effective strain rate below which the axes are
%                 treated as undefined and the result abstains (NaN).
%                 Principal directions from a near-zero strain rate are
%                 noise, and reporting them as a measurement is how a
%                 quiet divide turns into a spurious alignment.
%
% out.angle_deg  angle between the fabric axis and the reference principal
%                axis, in [0, 90], NaN where abstained
% out.reference  which principal axis was used
% out.ref_az     the reference azimuth actually used [rad]
% out.regime     passed through from strain when available
%
% See also ptt.strainRateAxes, ptt.orientationToFlow.
if nargin < 3, opts = struct(); end
tf = H_req(opts, 'theta_frame', ...
  'opts.theta_frame must be ''true'' or ''grid''');
ref_kind = lower(char(H_opt(opts, 'reference', 'extension')));
if ~ismember(ref_kind, {'extension', 'compression'})
  error('ptt:orientationToStrain:reference', ...
    'opts.reference must be ''extension'' or ''compression''');
end

regime = [];
e_eff = [];
if isstruct(strain)
  if ~isfield(strain, 'az_ext')
    error('ptt:orientationToStrain:strain', ...
      'strain struct must come from ptt.strainRateAxes (needs .az_ext)');
  end
  if strcmp(ref_kind, 'extension'), ref_az = strain.az_ext; else, ref_az = strain.az_cmp; end
  if isfield(strain, 'regime'), regime = strain.regime; end
  if isfield(strain, 'e_eff'), e_eff = strain.e_eff; end
  sf = lower(char(H_opt(opts, 'strain_frame', 'grid')));
else
  ref_az = double(strain);
  sf = lower(char(H_opt(opts, 'strain_frame', 'grid')));
end
if ~ismember(sf, {'true', 'grid'})
  error('ptt:orientationToStrain:frame', ...
    'opts.strain_frame must be ''true'' or ''grid''');
end

if ~strcmp(tf, sf)
  if ~isfield(opts, 'convergence') || isempty(opts.convergence)
    error('ptt:orientationToStrain:convergence', ...
      ['theta0 is in the %s frame and the strain axes are in the %s frame; ' ...
       'supply opts.convergence (grid-to-true, rad) rather than differencing them'], ...
      tf, sf);
  end
  cvg = double(opts.convergence);
  if strcmp(sf, 'grid'), ref_az = ref_az + cvg; else, ref_az = ref_az - cvg; end
end

ang = axisFold(theta0, ref_az);

e_min = H_opt(opts, 'e_min', 0);
if ~isempty(e_eff) && e_min > 0
  ang(e_eff < e_min) = NaN;
end

out = struct();
out.angle_deg = ang;
out.reference = ref_kind;
out.ref_az = ref_az;
out.regime = regime;
end

function v = H_req(s, f, msg)
if ~isstruct(s) || ~isfield(s, f) || isempty(s.(f))
  error('ptt:orientationToStrain:frame', '%s', msg);
end
v = lower(char(s.(f)));
if ~ismember(v, {'true', 'grid'})
  error('ptt:orientationToStrain:frame', '%s', msg);
end
end

function v = H_opt(s, f, d)
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
