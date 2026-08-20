%TEST_QUADPOL_SEGMENTED The laterally-segmented first pass: superset + rescue.
%
% Validates ptt.quadpolFrameTheta + ptt.thetaProfileAt - the segmented
% frame pass the pipeline hands to its per-block section fits - against
% the two properties the design must have:
%
%   A SUPERSET (gate a): on a laterally-uniform frame the segmented pass
%     must reproduce the current architecture. Old handoff (one pooled
%     frame theta0 profile) and new handoff (per-segment profiles,
%     interpolated per block) must give the same block dlam to 0.01, and
%     the segment axes must track truth.
%
%   THE RESCUE (gate b): test_egrip_blocks.m proved the killer: a 0.10
%     along-frame dlam ramp kills the pooled frame pass by ~100 m depth
%     and NO block size rescues the section, because every block inherits
%     the dead frame theta0. The segmented pass must rescue exactly this
%     case: block dlam med|err| < 0.03 against the local truth and
%     correlation > 0.7 along the ramp.
%
%   LATERAL AXIS (the margin case): a mid-frame 45 deg axis step (what a
%     shear margin does, and what the EastGRIP movie shows) must appear
%     in the segment profiles - outer segments within 3 deg of their
%     truths - while block dlam stays unbiased.
%
% Environment: the test_egrip_blocks forward (flat SNR, standard
% correlated + decorrelated pedestal, 9 m traces, 14-trace blocks, cap
% 0.45) so case B is the documented failure, bit-comparable.
%
% Run: matlab -batch "run('opr_fabric/test/test_quadpol_segmented.m')"
clear;
rng(7);
t0 = tic;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', '..'));

C = ptt.constants();
fc = 750e6;
gpd = 2*pi*fc*C.deps / (sqrt(C.eps_bar) * C.c*1e9);

Nz = 2800; dz = 0.5;
z = (0:Nz-1).' * dz;
Nx = 504; dxs = 9.0;
x = (0:Nx-1) * dxs;
NBLK = 14;
LEAK_C = 0.30 * exp(0.7i);
LEAK_D = 0.35;
NA = 0.5;

% helper options: test-speed grids; segment length 1.5 km -> 3 segments
FOPTS = struct('fc', fc, 'deramped', false, 'dlam_max', 0.45, ...
  'psi_step_deg', 4, 'theta_step_deg', 4, 'step_m', 30, ...
  'psi_step_seg_deg', 4, 'seg_len_m', 1500, 'track_az', 0);
% block-fit options (theta0 handed in per path below)
BOPTS = struct('fc', fc, 'deramped', false, 'dlam_max', 0.45, ...
  'psi_step_deg', 4, 'theta_step_deg', 4, 'step_m', 30);
az_tr = zeros(1, Nx);          % straight line, antennas = geographic

cases = { ...
  'A uniform (superset)', 0.32 + zeros(1, Nx), deg2rad(35) + zeros(1, Nx); ...
  'B dlam ramp (rescue)', 0.25 + 0.10 * (x / x(end)), deg2rad(35) + zeros(1, Nx); ...
  'C axis step (margin)', 0.32 + zeros(1, Nx), deg2rad(35 + 22.5) - deg2rad(45) * (x < x(end)/2)};

fails = 0;
for ci = 1:size(cases, 1)
  [nm, dlx, thx] = cases{ci, :};
  rng(40 + ci);
  r = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  g = (randn(Nz, Nx) + 1i*randn(Nz, Nx)) / sqrt(2);
  S = struct('hh', zeros(Nz, Nx), 'vv', zeros(Nz, Nx), ...
    'hv', zeros(Nz, Nx), 'vh', zeros(Nz, Nx));
  for j = 1:Nx
    ex = exp(1i * gpd * dlx(j) * z);
    cd_ = cos(thx(j)); sd = sin(thx(j));
    hh = cd_^2 * ex + sd^2;
    vv = sd^2 * ex + cd_^2;
    hv = cd_*sd * (ex - 1);
    S.hh(:, j) = hh .* r(:, j);
    S.vv(:, j) = vv .* r(:, j);
    S.hv(:, j) = hv .* r(:, j) + LEAK_C * ((hh + vv)/2) .* r(:, j) ...
      + LEAK_D * g(:, j);
    S.vh(:, j) = S.hv(:, j);
  end
  for f = {'hh','vv','hv','vh'}
    S.(f{1}) = S.(f{1}) + NA*(randn(Nz,Nx)+1i*randn(Nz,Nx));
  end

  % --- the segmented frame pass
  fp = ptt.quadpolFrameTheta(S, z, az_tr, x, FOPTS);

  % segment-axis recovery against the segment-mean truth
  seg_err = nan(1, fp.nseg);
  for s = 1:fp.nseg
    js = x >= fp.seg_x(s) - 750 & x < fp.seg_x(s) + 750;
    th_true_s = angle(mean(exp(2i*thx(js)))) / 2;
    m = fp.zw > 300 & fp.zw < 1100 & isfinite(fp.th_seg(:, s));
    if ~any(m), continue; end
    dth = mod(rad2deg(fp.th_seg(m, s) - th_true_s) + 90, 180) - 90;
    seg_err(s) = abs(median(dth));
  end

  % --- block pass twice: OLD handoff (pooled frame profile) and NEW
  okf = isfinite(fp.lsq.theta0);
  th_prof_old = struct('z', fp.lsq.zw(okf), 'theta', ...
    fp.lsq.theta0(okf) + deg2rad(FOPTS.track_az));
  % the OLD path hands blocks the UNSMOOTHED pooled frame theta0 plus the
  % pooled frame pedestal - the test_egrip_blocks.m killer convention, so
  % gate b reproduces the documented failure (the q-weighted smoothing of
  % fp.th_frame bridges this synthetic ramp and would hide it); the NEW
  % path takes both from the segmented pass (fp.ped_ant may be
  % segment-derived)
  ped_old = fp.lsq.pedestal;
  ped_new = fp.ped_ant;
  nblk = floor(Nx / NBLK);
  dl_old = nan(1, nblk); dl_new = nan(1, nblk); tr_blk = nan(1, nblk);
  for b = 1:nblk
    j0 = (b-1)*NBLK + 1; j1 = min(b*NBLK, Nx);
    Sb = struct('hh', S.hh(:, j0:j1), 'vv', S.vv(:, j0:j1), ...
      'hv', S.hv(:, j0:j1), 'vh', S.vh(:, j0:j1));
    xb = mean(x(j0:j1));
    o = BOPTS;
    if all(isfinite(ped_old)), o.pedestal = ped_old; end
    o.theta0 = th_prof_old;
    ob = ptt.quadpolFabricLS(Sb, z, o);
    m = ob.zw > 300 & ob.zw < 1100;
    dl_old(b) = median(ob.dlam(m), 'omitnan');
    o = BOPTS;
    if all(isfinite(ped_new)), o.pedestal = ped_new; end
    o.theta0 = ptt.thetaProfileAt(fp, xb);   % geographic; heading 0
    ob = ptt.quadpolFabricLS(Sb, z, o);
    dl_new(b) = median(ob.dlam(m), 'omitnan');
    tr_blk(b) = mean(dlx(j0:j1));
  end
  ok = isfinite(dl_new) & isfinite(tr_blk);
  err_new = median(abs(dl_new(ok) - tr_blk(ok)));
  oko = isfinite(dl_old) & isfinite(tr_blk);
  err_old = median(abs(dl_old(oko) - tr_blk(oko)));
  if nnz(ok) >= 3 && std(tr_blk(ok)) > 1e-6
    cc = corrcoef(dl_new(ok), tr_blk(ok)); corr_new = cc(1, 2);
  else
    corr_new = NaN;
  end
  dl_delta = median(abs(dl_new(ok & oko) - dl_old(ok & oko)));

  fprintf(['%s: nseg %d | seg th err [%s] deg | blk dlam err old %.3f ' ...
    'new %.3f | corr %.2f | old-new %.4f\n'], nm, fp.nseg, ...
    strtrim(sprintf('%.1f ', seg_err)), err_old, err_new, corr_new, dl_delta);

  switch ci
    case 1   % superset: new == old == truth
      okc = dl_delta < 0.01 && err_new < 0.02 && max(seg_err) < 1.5;
    case 2   % rescue: the documented killer must be fixed by segmentation
      okc = err_new < 0.03 && corr_new > 0.7 && err_old > 0.05;
    case 3   % margin: outer segments hit their axes, dlam unbiased
      okc = seg_err(1) < 3 && seg_err(end) < 3 && err_new < 0.03;
  end
  fprintf('  -> %s\n', H_tick(okc));
  fails = fails + ~okc;
end

fprintf('\n%s (%.1f min)\n', H_tick(fails == 0), toc(t0)/60);
if fails > 0
  error('test_quadpol_segmented:failed', '%d case(s) failed', fails);
end

function s = H_tick(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
