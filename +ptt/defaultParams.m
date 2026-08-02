function par = defaultParams()
%DEFAULTPARAMS Default firn-ice column parameters (Rathmann 2026, fig. 5/6).
%
%   Coordinate convention: z is HEIGHT above the bed, so zhat = z/H is 0 at
%   the bed and 1 at the surface.
%
%   Fields:
%     H         total ice-column thickness [m]
%     rho_sfc   relative density of surface snow (Herron-Langway, eq. 3.5)
%     rho_bco   relative density at bubble close-off
%     zhat_bco  relative height z/H of bubble close-off
%     e0        near-surface bubble eccentricity (eq. 3.4)
%     p         eccentricity power-law exponent
%     lam_x_sfc, lam_z_sfc   fabric eigenvalues at the surface
%     lam_x_bed, lam_z_bed   fabric eigenvalues at the bed
%                            (lam_y = 1 - lam_x - lam_z, eq. 4.2)
%     nz        number of grid points along a ray path

par.H        = 250;
par.rho_sfc  = 0.35;
par.rho_bco  = 0.81;
par.zhat_bco = 0.9;
par.e0       = 0.8;
par.p        = 0.5;

par.lam_x_sfc = 1/3;
par.lam_z_sfc = 1/3;
par.lam_x_bed = 1/3;
par.lam_z_bed = 1/3;

par.nz = 2001;

end
