function [dtau, theta0] = twttDifference(par, L, z)
%TWTTDIFFERENCE Oblique two-way traveltime difference between polarizations.
%   [dtau, theta0] = TWTTDIFFERENCE(par, L, z) models the TWTT difference
%   dtau = t_xz - t_y [ns] between xz- and y-polarized radar waves for a
%   common-midpoint geometry with half-offset L [m] reflecting at height z
%   [m] above the bed (Rathmann 2026, eqs. 2.11-2.13 and appendix A).
%   Refractive bending in firn is included by solving for the ray parameter
%   that connects transmitter and midpoint, then integrating slowness along
%   the bent path. theta0 is the initial propagation angle [deg].
%
%   L and z may be arrays of equal size. Geometries inside the refractive
%   shadow zone (unreachable for any launch angle) return NaN.
%
%   Notes:
%   - The x-axis is the horizontal fabric principal axis in the survey
%     plane; antennas are assumed aligned with the fabric principal frame.
%   - Both polarizations are propagated along the y-wave ray path; solid
%     ice is so weakly birefringent that the paths are virtually identical.

sz = size(L);
L = L(:); z = z(:);
assert(numel(L) == numel(z), 'L and z must have the same number of elements.');
assert(all(z >= 0 & z < par.H), 'Reflector heights must satisfy 0 <= z < H.');

n      = numel(L);
dtau   = nan(n, 1);
theta0 = nan(n, 1);

for i = 1:n
    zpath = linspace(z(i), par.H, par.nz).';
    P  = ptt.columnProfiles(par, zpath / par.H);
    Sy = P.S(:, 2);
    Vx = P.V(:, 1);
    Vz = P.V(:, 3);
    S0 = Sy(end); % surface slowness (slowness decreases with height)

    if L(i) <= 0
        q = 0;
    else
        % Ray parameter q = sin(theta0)*S0 from the offset integral (A1),
        % monotonically increasing in q; solve by bisection. u -> 1 at the
        % surface as q -> min(Sy), where the offset integral diverges.
        offsetFun = @(q) trapz(zpath, offsetIntegrand(q, Sy));
        qhi = (1 - 1e-12) * min(Sy);
        if offsetFun(qhi) < L(i)
            continue % refractive shadow zone: geometry unreachable
        end
        qlo = 0;
        for it = 1:80
            qm = (qlo + qhi) / 2;
            if offsetFun(qm) < L(i)
                qlo = qm;
            else
                qhi = qm;
            end
        end
        q = (qlo + qhi) / 2;
    end

    % Snell angle along the path (eq. 2.13) and slowness of the
    % xz-polarized wave (eq. 2.12); geometric factor g = 1/cos(theta).
    u    = min(q ./ Sy, 1 - 1e-15);
    cos2 = 1 - u.^2;
    cth  = sqrt(cos2);
    Sxz  = 1 ./ sqrt(cos2 .* Vx.^2 + u.^2 .* Vz.^2);

    dtau(i)   = 2 * trapz(zpath, (Sxz - Sy) ./ cth);
    theta0(i) = asind(q / S0);
end

dtau   = reshape(dtau, sz);
theta0 = reshape(theta0, sz);

end

function f = offsetIntegrand(q, Sy)
u = min(q ./ Sy, 1 - 1e-12);
f = u ./ sqrt(1 - u.^2); % tan(theta(z)), eq. A1
end
