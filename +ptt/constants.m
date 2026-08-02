function C = constants()
%CONSTANTS Physical constants for the polarimetric traveltime model.
%   Units: meters and nanoseconds throughout, permittivities relative.
%   Ice values from Rathmann (2026) eq. (2.3), valid for ice-sheet
%   temperatures and radar frequencies (Fujita et al., 2000).

C.c       = 0.299792458; % vacuum speed of light [m/ns]
C.eps_bar = 3.171;       % average relative permittivity of an ice grain
C.deps    = 0.034;       % grain permittivity anisotropy eps_c - eps_a
C.eps_air = 1.0;         % relative permittivity of air

end
