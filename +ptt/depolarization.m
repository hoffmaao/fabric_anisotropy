function Nz = depolarization(e)
%DEPOLARIZATION Vertical depolarization factor of a prolate spheroid.
%   Nz = DEPOLARIZATION(e) evaluates Rathmann (2026) eq. (3.1) for standing
%   (vertically elongated) bubbles with eccentricity e in [0, 1).
%   The horizontal factors follow from the unit trace: Nx = Ny = (1-Nz)/2.

Nz = zeros(size(e));

% Series expansion near e = 0 avoids 0/0: Nz = 1/3 - 2e^2/15 + O(e^4)
small = e < 1e-3;
Nz(small) = 1/3 - 2/15*e(small).^2;

es = e(~small);
Nz(~small) = (1 - es.^2) ./ (2*es.^3) .* (log((1 + es)./(1 - es)) - 2*es);

end
