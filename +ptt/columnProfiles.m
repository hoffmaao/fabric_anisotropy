function P = columnProfiles(par, zhat)
%COLUMNPROFILES Dielectric profiles of the firn-ice column.
%   P = COLUMNPROFILES(par, zhat) evaluates, at relative heights zhat = z/H
%   (0 = bed, 1 = surface), the closed-form profiles of Rathmann (2026):
%     - relative density rho: single-stage Herron-Langway model (eq. 3.5)
%     - bubble eccentricity e: power law vanishing at 0.99*zhat_bco (eq. 3.4)
%     - fabric eigenvalues lam = [lam_x lam_y lam_z]: linear in z (eq. 4.2)
%     - solid-ice eigenpermittivities (eq. 2.5)
%     - firn eigenpermittivities eps via Maxwell-Garnett mixing (eq. 3.2)
%       with the bubble depolarization tensor (eq. 3.1)
%     - eigenvelocities V [m/ns] and eigenslownesses S [ns/m] (eq. 2.6)
%   Columns of lam, eps, V, S correspond to the (x, y, z) principal axes.

C = ptt.constants();
zhat = zhat(:);

% Herron-Langway single-stage density (eq. 3.5)
beta  = log(1/par.rho_sfc - 1);
alpha = (log(1/par.rho_bco - 1) - beta) / (1 - par.zhat_bco);
rho   = 1 ./ (1 + exp(alpha*(1 - zhat) + beta));

% Bubble eccentricity power law, vanishing just below close-off (eq. 3.4)
zhat_fb = 0.99 * par.zhat_bco;
zeta    = max((zhat - zhat_fb) / (1 - zhat_fb), 0);
e       = par.e0 * zeta.^par.p;

% Fabric eigenvalue profiles, lam_y by normalization. Default: linear
% between surface and bed values (eq. 4.2). If par.fabric is set, it
% overrides with node-based profiles: fields zhat (ascending), lam_x,
% lam_z, and optional interpolation method (default 'linear'; use
% 'previous' for piecewise-constant intervals with nodes at the interval
% bottom edges). Queries are clamped to the node range.
if isfield(par, 'fabric') && ~isempty(par.fabric)
    f = par.fabric;
    if ~isfield(f, 'method'), f.method = 'linear'; end
    zq = min(max(zhat, min(f.zhat)), max(f.zhat));
    lam_x = interp1(f.zhat(:), f.lam_x(:), zq, f.method);
    lam_z = interp1(f.zhat(:), f.lam_z(:), zq, f.method);
else
    lam_x = par.lam_x_bed + (par.lam_x_sfc - par.lam_x_bed) .* zhat;
    lam_z = par.lam_z_bed + (par.lam_z_sfc - par.lam_z_bed) .* zhat;
end
lam = [lam_x, 1 - lam_x - lam_z, lam_z];

% Solid-ice eigenpermittivities (eq. 2.5)
eps_ice = C.eps_bar + C.deps * (lam - 1/3);

% Maxwell-Garnett mixing with air inclusions (eqs. 3.1, 3.2)
Nz  = ptt.depolarization(e);
N   = [(1 - Nz)/2, (1 - Nz)/2, Nz];
phi = (1 - rho) ./ (1 + N .* rho .* (C.eps_air ./ eps_ice - 1));
eps = eps_ice + phi .* (C.eps_air - eps_ice);

P.zhat = zhat;
P.rho  = rho;
P.e    = e;
P.lam  = lam;
P.eps  = eps;
P.V    = C.c ./ sqrt(eps);
P.S    = 1 ./ P.V;

end
