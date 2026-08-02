function [dlam, out] = invertHorizontalFabric(obs, par)
%INVERTHORIZONTALFABRIC Horizontal fabric contrast from common-offset data.
%   [dlam, out] = INVERTHORIZONTALFABRIC(obs, par) inverts polarimetric
%   TWTT differences collected at a fixed (common) offset for the profile
%   of the horizontal fabric contrast
%       dlam(z) = lam_x(z) - lam_y(z),
%   the component of the orientation tensor that common-offset (nadir-type)
%   sounding is sensitive to (Rathmann 2026, eq. 2.10). dlam is assumed
%   piecewise constant between consecutive picked reflectors and solved
%   top-down by layer stripping: each layer's TWTT difference increment
%   determines the contrast of the interval directly above it via a
%   bisection solve through the full Maxwell-Garnett forward model, so the
%   density and bubble corrections and (small) offset obliquity are all
%   accounted for.
%
%   obs is a struct with fields:
%       L     common half-offset [m] (scalar; 0 for nadir)
%       z     reflector heights above the bed [m], any order
%       dtau  cumulative TWTT differences t_x - t_y [ns] from the surface
%             down to each reflector (positive when lam_x > lam_y)
%   par is the model setup (see ptt.defaultParams): H, density and bubble
%   parameters, plus an ASSUMED vertical eigenvalue profile lam_z(z) given
%   by lam_z_sfc/lam_z_bed (common-offset data cannot constrain lam_z; use
%   1/3 near the surface or values from a nearby core/CMP survey).
%
%   Outputs:
%       dlam        inferred contrast per interval, ordered like out.ztop
%       out.ztop    interval top heights [m] (first = surface)
%       out.zbot    interval bottom heights [m] (the reflectors)
%       out.dtau_fit modeled TWTT differences at the reflectors [ns]
%       out.par     par with the inferred fabric attached (par.fabric)
%
%   The split between lam_x and lam_y uses the horizontal budget
%   h(z) = 1 - lam_z(z): lam_x = (h + dlam)/2, lam_y = (h - dlam)/2,
%   which bounds |dlam| <= h.

[zs, order] = sort(obs.z(:), 'descend'); % shallowest reflector first
dts = obs.dtau(:);
dts = dts(order);
m   = numel(zs);
assert(all(zs > 0 & zs < par.H), 'Reflector heights must satisfy 0 < z < H.');

zhat_nodes = [zs(end:-1:1) / par.H; 1]; % ascending: deepest ... shallowest, surface

% Assumed lam_z, evaluated at interval midpoints (step approximation)
zhat_mid = [(1 + zs(1)/par.H)/2; (zs(1:end-1) + zs(2:end)) / (2*par.H)];
lam_z_mid = par.lam_z_bed + (par.lam_z_sfc - par.lam_z_bed) .* zhat_mid;
h_mid = 1 - lam_z_mid; % horizontal eigenvalue budget per interval

    function par2 = buildPar(dl)
        % Piecewise-constant fabric: node at each interval's bottom edge
        % ('previous' interpolation + clamping covers up to the surface).
        % dl(k) is the contrast of interval k (k = 1 shallowest).
        k = numel(dl);
        dl_nodes = [repmat(dl(k), m - k + 1, 1); dl(k-1:-1:1); dl(1)];
        hz = [repmat(h_mid(k), m - k + 1, 1); h_mid(k-1:-1:1); h_mid(1)];
        lz = [repmat(lam_z_mid(k), m - k + 1, 1); lam_z_mid(k-1:-1:1); lam_z_mid(1)];
        par2 = par;
        par2.fabric = struct('method', 'previous', 'zhat', zhat_nodes, ...
            'lam_x', (hz + dl_nodes)/2, 'lam_z', lz);
    end

dlam = zeros(m, 1);
dtau_fit = zeros(m, 1);
for k = 1:m
    % Solve F(dlam_k) = dts(k) by bisection; F is monotone increasing.
    lo = -h_mid(k) + 1e-6;
    hi =  h_mid(k) - 1e-6;
    Ffun = @(dl_k) ptt.twttDifference(buildPar([dlam(1:k-1); dl_k]), obs.L, zs(k));
    Flo = Ffun(lo); Fhi = Ffun(hi);
    if ~(isfinite(Flo) && isfinite(Fhi))
        error('ptt:shadowZone', 'Reflector %d unreachable for L = %g m.', k, obs.L);
    end
    if dts(k) <= Flo
        dlam(k) = lo;
    elseif dts(k) >= Fhi
        dlam(k) = hi;
    else
        for it = 1:60
            mid = (lo + hi) / 2;
            if Ffun(mid) < dts(k), lo = mid; else, hi = mid; end
        end
        dlam(k) = (lo + hi) / 2;
    end
    dtau_fit(k) = Ffun(dlam(k));
end

out.ztop     = [par.H; zs(1:end-1)];
out.zbot     = zs;
out.dtau_fit = dtau_fit;
out.par      = buildPar(dlam);

end
