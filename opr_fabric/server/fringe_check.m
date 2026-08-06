%FRINGE_CHECK Is dlam proportional to the observed fringe density?
%
% Fringe density per unit TWTT is fc*d(dtau)/dt, and the forward model
% makes d(dtau)/dt a function of dlam alone (through the birefringence of
% the firn-ice column). So the interferogram gives dlam directly, with no
% inversion involved, and the inverted profile has to reproduce it. This
% checks that end to end.
%
% Every constant is taken from ptt itself rather than hand-derived: the
% conversion is measured by perturbing dlam in ptt.columnProfiles and
% reading the change in the y/x eigen-slowness difference, so whatever
% Rathmann's permittivities and the Maxwell-Garnett firn mixing actually
% do is what is used here.

scratch = '/kucresis/scratch/hoffmana_sta/fabric';
in_fn = ['/cresis/nvme/opr_data/accum/2024_Antarctica_Ground2/' ...
  'CSARP_polarimetric/20250108_02/Data_20250108_02_009.mat'];
sec_fn = fullfile(scratch, 'stages', 'fabric_sections.mat');
out_fn = fullfile(scratch, 'stages', 'fringe_check.mat');

code = fullfile(scratch, 'code');
addpath(code); addpath(fullfile(code,'opr_fabric'));
addpath(fullfile(code,'opr_fabric','test','stubs'));

FC = 750e6;
par = ptt.defaultParams();
par.H = 2000; par.lam_z_sfc = 1/3; par.lam_z_bed = 1/3;
par.zhat_bco = 1 - 60/par.H;

%% 1. ptt's own conversion: d(dtau)/dt per unit dlam, versus depth
% ptt.columnProfiles returns slowness S in ns/m for the x and y eigen-
% directions. d(dtau)/dz = (S_y - S_x) [ns/m] and dz/dt = 1/(2*S_mean),
% so d(dtau)/dt = (S_y - S_x) / (2*S_mean) - dimensionless, and the
% fringe rate is fc times that.
% Call the forward model the inversion itself uses, with a uniform dlam,
% rather than deriving the constant by hand - the first attempt at this
% dropped the two-way factor that lives inside ptt.twttDifference
% (dtau = 2*trapz(...)) and came out 2x high, which is exactly the size
% of discrepancy being investigated.
zhat = linspace(1, 0, 4001).';
depth = par.H * (1 - zhat);
dl_probe = 0.10;
p0 = par;
P0 = ptt.columnProfiles(p0, zhat);
twtt_col = 2e-9 * cumtrapz(depth, P0.S(:,2));

zq = (50:50:1950).';                 % depths to probe, m below surface
pf = par;
pf.fabric = struct('method', 'previous', ...
  'zhat', [0; 1], ...
  'lam_x', ((1 - 1/3) + dl_probe)/2 * [1; 1], ...
  'lam_z', (1/3) * [1; 1]);
% ptt.twttDifference returns dtau in NANOSECONDS (its slownesses are ns/m)
dtau_probe = 1e-9 * ptt.twttDifference(pf, zeros(size(zq)), par.H - zq);
% d(dtau)/d(twtt), dimensionless, per unit dlam
tq = interp1(depth, twtt_col, zq);
ddtau_dt_q = gradient(dtau_probe(:), tq(:)) / dl_probe;
ddtau_dt = interp1(zq, ddtau_dt_q, depth, 'linear', 'extrap');
fprintf('ptt conversion (fringe rate per unit dlam, cycles/us):\n');
for zprint = [100 300 600 1000 1500 1900]
  [~, k] = min(abs(depth - zprint));
  fprintf('  %5.0f m: %.3f  (d(dtau)/dt = %.6f per unit dlam)\n', ...
    zprint, FC*ddtau_dt(k)*1e-6, ddtau_dt(k));
end

%% 2. measured d(dtau)/dt from the validated SNAPHU phase
want = {'interferogram_mlook','interferogram_coherence','snaphu_out_phase', ...
  'Time','Surface','Latitude','Longitude'};
have = whos('-file', in_fn);
sel = intersect(want, {have.name});
pol = load(in_fn, sel{:});
t = pol.Time(:);
surf_med = median(pol.Surface(isfinite(pol.Surface)));
tb = t - surf_med;
coh = abs(pol.interferogram_coherence);

% coherence-weighted mean phase per fast-time bin, then a local slope
ph = pol.snaphu_out_phase;
wgt = double(coh);
phm = sum(ph .* wgt, 2) ./ max(sum(wgt, 2), eps);
cohm = median(coh, 2);

% independent, unwrap-free rate straight off the interferogram
D = 30;
X = pol.interferogram_mlook(1+D:end, :) .* conj(pol.interferogram_mlook(1:end-D, :));

%% 3. compare against the inverted profile
Sload = load(sec_fn); Sd = Sload.S.ridge_a;
zc = median((double(Sd.top) + double(Sd.bot))/2, 2, 'omitnan');
dl_inv = median(double(Sd.dlam), 2, 'omitnan');
q_inv = median(double(Sd.quality), 2, 'omitnan');

fprintf('\n%8s %8s %10s %10s %10s %8s\n', 'depth_m', 'twtt_us', ...
  'dlam_inv', 'dlam_snaphu', 'dlam_fringe', 'coh');
res = zeros(numel(zc), 5);
for i = 1:numel(zc)
  % twtt window spanning this interval
  z0 = median(double(Sd.top(i,:)), 'omitnan');
  z1 = median(double(Sd.bot(i,:)), 'omitnan');
  [~, k0] = min(abs(depth - z0));
  [~, k1] = min(abs(depth - z1));
  % twtt of those depths, from the same column model
  twtt = twtt_col;
  t0 = twtt(k0); t1 = twtt(k1);
  m = tb >= t0 & tb < t1;
  if nnz(m) < 8, continue; end
  % conversion at this interval's mid-depth
  [~, km] = min(abs(depth - (z0+z1)/2));
  conv = ddtau_dt(km);

  c = polyfit(t(m), phm(m), 1);
  rate_sn = c(1) / (2*pi);                    % cycles per second of twtt
  mm = m(1:end-D);
  r = sum(X(mm, :), 'all');
  rate_fr = angle(r) / (2*pi*D*(t(2)-t(1)));

  res(i,:) = [ (z0+z1)/2, 0.5*(t0+t1)*1e6, dl_inv(i), ...
    abs(rate_sn)/(FC*conv), abs(rate_fr)/(FC*conv) ];
  fprintf('%8.0f %8.2f %10.4f %10.4f %10.4f %8.2f\n', res(i,1), res(i,2), ...
    dl_inv(i), res(i,4), res(i,5), median(cohm(m)));
end

save(out_fn, '-v7', 'res', 'ddtau_dt', 'depth', 'zc', 'dl_inv', 'q_inv');

%% 4. profile azimuth - the eigenvalues are relative to THIS pair of axes
la = pol.Latitude(:); lo = pol.Longitude(:);
ok = isfinite(la) & isfinite(lo);
la = la(ok); lo = lo(ok);
dy = la(end) - la(1);
dx = (lo(end) - lo(1)) * cosd(mean(la));
az = mod(atan2d(dx, dy), 360);
fprintf('\nProfile heading %.1f deg E of N (start %.4f,%.4f -> end %.4f,%.4f)\n', ...
  az, la(1), lo(1), la(end), lo(end));
fprintf('  lam_along is the horizontal eigenvalue along that heading;\n');
fprintf('  lam_cross is along %.1f deg E of N (heading + 90).\n', mod(az+90,360));
fprintf('wrote %s\n', out_fn);
