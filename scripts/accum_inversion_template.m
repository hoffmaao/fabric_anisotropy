%ACCUM_INVERSION_TEMPLATE Horizontal fabric from accumulation-radar picks.
%   Template for the common-offset workflow: the CReSIS ground-based
%   accumulation radar is towed at a fixed Tx-Rx separation, and the two
%   co-polarized channels (or repeat passes with rotated antennas) give the
%   TWTT difference t_x - t_y of tracked reflectors. Layer stripping then
%   yields the horizontal fabric contrast dlam(z) = lam_x - lam_y per
%   interval (ptt.invertHorizontalFabric).
%
%   Conventions and assumptions:
%   - x is horizontal along the polarization/survey plane, y perpendicular;
%     these must be aligned with the horizontal fabric principal axes (find
%     the principal azimuth first, e.g. from azimuthal rotation tests, and
%     watch for polarization rotation if badly misaligned).
%   - dtau values are CUMULATIVE from the surface down to each reflector,
%     in ns, positive when the x-polarized wave is slower (lam_x > lam_y).
%   - Common-offset data cannot constrain lam_z or the BCO depth: the
%     density/bubble model and a lam_z profile are ASSUMED (defaults or a
%     nearby core). At small L the inferred dlam is insensitive to these
%     (see scripts/synthetic_common_offset.m). Extending to varying-offset
%     CMP acquisitions later unlocks the full inversion (ptt.invertFabric).

clear;
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..'));

%% ---- 1. Load picks (EDIT THIS) -----------------------------------------
H     = 2000;                % standing default; only needs to exceed the
                             % deepest reflector (twtt-depth map depends
                             % only on the close-off depth at nadir)
obs.L = 0;                   % half Tx-Rx separation; 0 for the ground accum
                             % radar (switch-gated crossed bowties, single
                             % phase center)
depth = [25; 50; 88; 125];   % reflector depths below surface [m]
obs.z = H - depth;           % heights above the bed
obs.dtau = nan(size(obs.z)); % <-- t_x - t_y [ns], cumulative per reflector

%% ---- 2. Assumed column model -------------------------------------------
par = ptt.defaultParams();
par.H = H;
par.zhat_bco = 1 - 60/H;     % bubble close-off at 60 m depth (site value)
% par.e0, par.p: keep defaults or set from a nearby core
% Assumed vertical eigenvalue profile (cannot be inferred from these data):
par.lam_z_sfc = 1/3;
par.lam_z_bed = 1/3;

%% ---- 3. Invert by layer stripping --------------------------------------
[dlam, out] = ptt.invertHorizontalFabric(obs, par);

fprintf('Interval (m a.b.)   dlam = lam_x - lam_y\n');
for k = 1:numel(dlam)
    fprintf('%6.0f - %5.0f      %8.3f\n', out.ztop(k), out.zbot(k), dlam(k));
end

%% ---- 4. Plot -----------------------------------------------------------
figure;
subplot(1, 2, 1);
stairs([dlam; dlam(end)], [out.zbot; H], 'LineWidth', 1.5);
xlabel('\Delta\lambda'); ylabel('z (m above bed)'); grid on;
subplot(1, 2, 2);
plot(obs.dtau, obs.z, 'ko', out.dtau_fit, obs.z, 'r^');
xlabel('\Delta\tau_{x,y} (ns)'); legend('observed', 'fit'); grid on;
