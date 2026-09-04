function [T, info] = coregisterChannels(S, ref, opts)
%COREGISTERCHANNELS Align polarimetric channels using the OPR coregistration.
%
% [T, info] = ptt.coregisterChannels(S, ref, opts)
%
% Thin wrapper around the CReSIS/OPR toolbox function `coregistration`
% (opr/matlab/processing/coregistration.m), applied to every channel of S
% against `ref` (default 'hh').
%
% CALLS THE TOOLBOX RATHER THAN REIMPLEMENTING IT. The production
% polarimetric products are built by polarimetric_task.m calling exactly
% this function, so using it here means the quad-pol path is aligned by the
% same code, with the same defaults, as the CSARP_polarimetric products
% every other result in this repo is derived from. A hand-rolled
% cross-spectrum fit was tried first and was worse on both counts: it
% reached |C_HHVV| 0.371 against the production 0.495, and it DEGRADED the
% already-aligned cross-polarized pair from 0.985 to 0.882, because a
% single global delay is the wrong model for an offset that varies over the
% image. The toolbox version estimates offsets on overlapping tiles in both
% range and azimuth, oversamples the correlation, median-filters the offset
% field and resamples by sinc interpolation.
%
% WHY IT MATTERS HERE. Reading the standardphase channels raw and skipping
% this step gives |C_HHVV| ~ 0.29 on Ridge A, against 0.495 for the
% coregistered product - and an earlier analysis in this repo drew a
% conclusion about antenna geometry from that uncorrected number. The
% channels are offset by a sub-sample instrument delay (row_offset ~ 0.78
% samples, 0.22 m), not by looking at different ice.
%
% Inputs
%   S     struct of complex [Nt x Nx] channels
%   ref   reference channel name (default 'hh')
%   opts  passed to `coregistration`; defaults below match polarimetric.m
%           Tt 51, Tx 101, overlap_t 25, overlap_x 25,
%           search_t 5, search_x 5, one_dim_search_en true
%         impl ('toolbox' default: the OPR function; 'ptt': the package
%           copy ptt.coregistration, identical unless its own options are
%           set) and, for 'ptt' only, peak_clamp_en (false) and
%           peak_contrast_min (0) - see ptt.coregistration for what they
%           fix and for the margin measurement that motivated the copy
%
% Output
%   T     struct of coregistered channels, ref passed through untouched
%   info.row_offset / .col_offset  per-channel offset fields
%   info.row_med                   median row offset per channel [samples]

if nargin < 2 || isempty(ref), ref = 'hh'; end
if nargin < 3, opts = struct(); end
impl = H_opt(opts, 'impl', 'toolbox');
if strcmpi(impl, 'ptt')
  coreg_fn = @ptt.coregistration;
  extra = {'peak_clamp_en', logical(H_opt(opts, 'peak_clamp_en', false)), ...
    'peak_contrast_min', H_opt(opts, 'peak_contrast_min', 0)};
else
  coreg_fn = @coregistration;
  extra = {};
end
if strcmpi(impl, 'toolbox') && exist('coregistration', 'file') ~= 2
  error('ptt:coregisterChannels:toolbox', ...
    ['the OPR toolbox function `coregistration` is not on the path; add ' ...
     'opr/matlab/processing (the startup script does this on the ' ...
     'server) - this wrapper deliberately does not carry its own ' ...
     'implementation']);
end
if ~isfield(S, ref)
  error('ptt:coregisterChannels:ref', 'reference channel %s not in S', ref);
end

% polarimetric.m's defaults, so this matches how the products were built
Tt = H_opt(opts, 'Tt', 51);
Tx = H_opt(opts, 'Tx', 101);
ov_t = H_opt(opts, 'overlap_t', 25);
ov_x = H_opt(opts, 'overlap_x', 25);
se_t = H_opt(opts, 'search_t', 5);
se_x = H_opt(opts, 'search_x', 5);
one_d = H_opt(opts, 'one_dim_search_en', true);

fn = fieldnames(S);
T = struct();
info = struct('row_offset', struct(), 'col_offset', struct(), ...
  'row_med', struct());
for k = 1:numel(fn)
  ch = fn{k};
  if strcmp(ch, ref)
    T.(ch) = S.(ch);
    info.row_med.(ch) = 0;
    continue;
  end
  [reg, row_off, col_off] = coregistration(S.(ref), S.(ch), ...
    'Tt', Tt, 'Tx', Tx, 'overlap_t', ov_t, 'overlap_x', ov_x, ...
    'search_t', se_t, 'search_x', se_x, 'one_dim_search_en', one_d, ...
    'debug_print_en', false);
  T.(ch) = reg;
  info.row_offset.(ch) = row_off;
  info.col_offset.(ch) = col_off;
  info.row_med.(ch) = median(row_off(isfinite(row_off)));
end

end

function v = H_opt(o, fld, d)
if isstruct(o) && isfield(o, fld) && ~isempty(o.(fld)), v = o.(fld);
else, v = d; end
end
