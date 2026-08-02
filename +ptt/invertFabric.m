function [parFit, out] = invertFabric(obs, parInit)
%INVERTFABRIC Infer fabric profile and bubble close-off depth from TWTT data.
%   [parFit, out] = INVERTFABRIC(obs, parInit) solves the inverse problem of
%   Rathmann (2026) sec. 4: find the five parameters
%       f = (zhat_bco, lam_x_sfc, lam_x_bed, lam_z_sfc, lam_z_bed)
%   minimizing the mean-square misfit J (eq. 4.1) between observed and
%   modeled TWTT differences over the sampled CMP geometries.
%
%   obs is a struct with equally sized fields:
%       L     half-offsets [m]
%       z     reflector heights above the bed [m]
%       dtau  observed TWTT differences t_xz - t_y [ns]
%   parInit provides both the fixed model setup (H, density, eccentricity)
%   and the initial guess for the five free parameters.
%
%   Uses fmincon (SQP, as in the paper) when the Optimization Toolbox is
%   available, otherwise falls back to fminsearch on logit-transformed
%   variables. Bounds keep all parameters in [0, 1]; linear constraints
%   keep lam_x + lam_z <= 1 (so lam_y >= 0) at both boundaries.

free = {'zhat_bco', 'lam_x_sfc', 'lam_x_bed', 'lam_z_sfc', 'lam_z_bed'};
x0   = cellfun(@(f) parInit.(f), free).';
dobs = obs.dtau(:);

    function par = assign(x)
        par = parInit;
        for k = 1:numel(free)
            par.(free{k}) = x(k);
        end
    end

    function J = cost(x)
        dmod = ptt.twttDifference(assign(x), obs.L(:), obs.z(:));
        r = dmod - dobs;
        r(isnan(r)) = 10; % penalize geometries the trial model cannot reach
        J = mean(r.^2);
    end

if exist('fmincon', 'file') == 2
    % lam_x + lam_z <= 1 at surface (x2 + x4) and bed (x3 + x5)
    A = [0 1 0 1 0;
         0 0 1 0 1];
    b = [1; 1];
    opts = optimoptions('fmincon', 'Algorithm', 'sqp', 'Display', 'off', ...
        'MaxFunctionEvaluations', 5000);
    [x, J, exitflag] = fmincon(@cost, x0, A, b, [], [], ...
        zeros(5, 1), ones(5, 1), [], opts);
else
    % Toolbox-free fallback: sigmoid keeps parameters in (0, 1), quadratic
    % penalty enforces lam_x + lam_z <= 1.
    sig = @(y) 1 ./ (1 + exp(-y));
    pen = @(x) 1e3 * (max(0, x(2) + x(4) - 1)^2 + max(0, x(3) + x(5) - 1)^2);
    y0 = log(x0 ./ (1 - x0));
    opts = optimset('MaxFunEvals', 1e4, 'MaxIter', 1e4, 'Display', 'off');
    [y, J, exitflag] = fminsearch(@(y) cost(sig(y)) + pen(sig(y)), y0, opts);
    x = sig(y);
end

parFit = assign(x);

out.x        = x(:).';
out.free     = free;
out.J        = J;
out.J0       = cost(x0);
out.exitflag = exitflag;

end
