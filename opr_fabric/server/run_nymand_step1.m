%RUN_NYMAND_STEP1 Nymand step 1 across a site, from the coregistration caches.
%
% Solves ONE fabric orientation and the scattering-ratio profile per frame
% from the azimuthal POWER anomalies, by the iterative linearisation of
% ptt.polarimetricInverse - the first step of the two-step scheme, under a
% constant eigenframe.
%
% WHY THIS IS WORTH RUNNING ON ITS OWN. It reaches theta0 by a route that
% shares no machinery with the production estimator: powers rather than
% the complex coherence, a Taylor-expanded least-squares solve rather than
% a grid search over a cost surface. Two independent methods on identical
% coregistered data is the strongest internal check available for the
% orientation question this project currently has open, where the Ridge A
% products aggregate to theta0_geo 151 deg against 96.5 deg on record.
%
% Loops inside ONE MATLAB session rather than launching per frame:
% the work is ~25 s a frame against a ~30 s MATLAB start, and parallel
% MATLAB launches on this node have hit licence-server failures.
%
% Variables expected in the caller (site_chain style):
%   fab_root   work root, e.g. /kucresis/scratch/hoffmana_sta/fabric
%   tags       cellstr of frame tags, or leave unset for FILTER below
%   filter     substring the cache name must contain (default '2025', the
%              Ridge A season)
% Optional: nr (range looks, 15), psi_step_deg (2), fc (750e6).
if ~exist('fab_root', 'var'), fab_root = '/kucresis/scratch/hoffmana_sta/fabric'; end
if ~exist('filter', 'var'), filter = '2025'; end
if ~exist('nr', 'var'), nr = 15; end
if ~exist('psi_step_deg', 'var'), psi_step_deg = 2; end
if ~exist('fc', 'var'), fc = 750e6; end

qdir = fullfile(fab_root, 'stages', 'quadpol');
cdir = fullfile(qdir, 'coreg_cache');
odir = fullfile(qdir, 'nymand');
if exist(odir, 'dir') ~= 7, mkdir(odir); end

if ~exist('tags', 'var') || isempty(tags)
  L = dir(fullfile(cdir, ['creg_*' filter '*.mat']));
  tags = cell(numel(L), 1);
  for i = 1:numel(L)
    n = L(i).name; tags{i} = n(6:end-4);
  end
  tags = sort(tags);
end
fprintf('%d frame(s) matching filter "%s"\n\n', numel(tags), filter);

psi = (0:psi_step_deg:(180-psi_step_deg)) * pi/180;
an = @(P) 20*log10(max(P, realmin) ./ max(mean(P, 2, 'omitnan'), realmin));

nok = 0;
for i = 1:numel(tags)
  tag = tags{i};
  ofn = fullfile(odir, ['nymand_step1_' tag '.mat']);
  if exist(ofn, 'file') == 2
    fprintf('[%2d/%2d] %s  already done\n', i, numel(tags), tag); nok = nok + 1; continue
  end
  cfn = fullfile(cdir, ['creg_' tag '.mat']);
  sfn = fullfile(qdir, ['quadpol_section_' tag '.mat']);
  if exist(cfn, 'file') ~= 2 || exist(sfn, 'file') ~= 2
    fprintf('[%2d/%2d] %s  SKIP (cache or section missing)\n', i, numel(tags), tag); continue
  end
  t0 = tic;
  try
    c = load(cfn);
    z = double(c.z(:));
    S = struct('hh', c.hh, 'vv', c.vv, 'hv', c.hv, 'vh', c.vh);
    M = ptt.quadpolMoments(S, [nr size(S.hh, 2)]);
    clear c S
    A = ptt.quadpolAzimuth(M, psi);
    clear M

    % Power anomalies in the convention ptt.polarimetricInverse expects:
    % 20*log10 of the power ratio to the azimuthal mean, matching the
    % generator in test_polarimetric_inverse.
    obs = struct('psi', psi, 'dP_hh', an(double(A.Phh)), ...
      'dP_hv', an(double(A.Pxc)));
    clear A

    % The contrast for step 1's phase term comes from the frame's existing
    % quadpol profile. Step 1 needs a dlam it does not yet have, and the
    % scheme's step 2 supplies it from travel-time anomalies - which needs
    % a delay measured on UNREGISTERED channels, so it cannot come from
    % these caches. Bootstrapping from the production profile keeps step 1
    % honest about the phase while that pass is built; theta0 is what this
    % run is for, and it is what step 1 determines.
    sec = load(sfn); r = sec.res;
    zs = double(r.z(:)); dl = double(r.dlam_ls(:));
    n = min(numel(zs), numel(dl)); zs = zs(1:n); dl = dl(1:n);
    gd = isfinite(dl);
    if nnz(gd) < 5, error('no finite dlam in the section product'); end
    dl_z = interp1(zs(gd), dl(gd), z, 'linear', 'extrap');
    dl_z(~isfinite(dl_z)) = median(dl(gd));

    o = ptt.polarimetricInverse(obs, z, struct('fc', fc, 'dlam', dl_z, ...
      'n_r', 10, 'eta', 1e-2, 'azimuth_source', 'synthetic'));

    % production answer on its own grid, same trusted band, circular median
    th_ls = double(r.theta0_ls(:));
    m = min(numel(zs), numel(th_ls));
    band = zs(1:m) >= 200 & zs(1:m) <= 1200;
    gl = band & isfinite(th_ls(1:m));
    if nnz(gl) >= 3
      th_prod = angle(median(exp(2i*th_ls(gl))))/2;
    else
      th_prod = NaN;
    end
    sep = rad2deg(abs(angle(exp(2i*(o.theta0 - th_prod)))/2));

    res = struct();
    res.tag = tag;
    res.theta0 = o.theta0;            % antenna frame, synthesis sense
    res.theta0_prod = th_prod;        % same frame and sense
    res.sep_deg = sep;
    res.r_z = o.r_z; res.z = z;
    res.n_iter = o.n_iter;
    res.dlam_in = dl_z;
    if isfield(r, 'track_az'), res.track_az = r.track_az; end
    if isfield(r, 'lat'), res.lat = r.lat; res.lon = r.lon; end
    if isfield(r, 'ls_pedestal'), res.ls_pedestal = r.ls_pedestal; end
    tmp = [ofn '.tmp'];
    save(tmp, '-v7.3', 'res'); movefile(tmp, ofn);
    nok = nok + 1;
    fprintf('[%2d/%2d] %s  nymand %6.1f  prod %6.1f  sep %5.1f deg  (%.0f s)\n', ...
      i, numel(tags), tag, mod(rad2deg(o.theta0),180), ...
      mod(rad2deg(th_prod),180), sep, toc(t0));
  catch err
    fprintf('[%2d/%2d] %s  FAIL %s (%.0f s)\n', i, numel(tags), tag, err.message, toc(t0));
  end
end
fprintf('\n%d of %d frame(s) produced a step-1 result\n', nok, numel(tags));
