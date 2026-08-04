function A_full = cellToGrid(A_cell, t_cell, x_cell, Time, Nx)
%CELLTOGRID Interpolate an analysis-cell field onto the full data grid.
%   A_full = CELLTOGRID(A_cell, t_cell, x_cell, Time, Nx) resamples an
%   (ntc x nxc) quantity defined on the delta-k analysis-cell grid onto the
%   (Nt x Nx) fast-time / along-track grid the rest of the chain works on,
%   by bilinear interpolation in cell-index space with edge clamping.
%
%   ptt.deltakTraveltime uses this for the dtau it returns; the per-stage
%   ladder diagnostic (opr_fabric/server/run_deltak_stages.m) uses it to
%   put stages A and Q on that same grid, so every rung can be reduced by
%   the identical ptt.blockAverage call rather than by a separate
%   cell-column statistic. Comparing a stage reduced one way against an
%   estimator reduced another would confound the reduction with the
%   amplitude difference the diagnostic is trying to measure.

ntc = size(A_cell,1);
nxc = size(A_cell,2);

ti = interp1(t_cell, (1:ntc).', Time(:), 'linear', 'extrap');
ti = min(max(ti, 1), ntc);
xi = interp1(x_cell, (1:nxc).', (1:Nx).', 'linear', 'extrap');
xi = min(max(xi, 1), nxc);

A_full = interp2(A_cell, xi.', ti, 'linear');

end
