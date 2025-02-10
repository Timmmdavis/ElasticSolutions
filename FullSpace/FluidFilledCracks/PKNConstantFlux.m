function [L,  w, p, dpdx, v_x] = PKNConstantFlux(Q0, t, Eprime, nu, etaprime, H, x)
% PKNConstantFlux: Computes fracture half-length, fracture width
% profile, and net pressure for the storage-dominated (M-regime)
% similarity solution based on Kovalyshen & Detournay (2010).
%
% Arguments (input):
%   Q0       - Constant volumetric flow rate (m^3/s).
%   t        - Time (s).
%   Eprime   - Plane strain Young's modulus (Pa).
%   nu      - Poisson's ratio (dimensionless)#
%   etaprime - Effective viscosity parameter (Pa·s).
%   H        - Fracture height (m).
%   N        - Number of spatial points for discretization.
%
% Arguments (output):
%   L    - Fracture half-length (m).
%   w    - Fracture width profile (m).
%   p    - Net pressure along the fracture (Pa).
%   dpdx - Pressure gradient along the fracture (Pa/m).
%   v_x  - Crack tip velocity (m/s).
%
% References:
% Kovalyshen & Detournay (2010) - https://doi.org/10.1007/s11242-009-9403-4
%
% Example usage:
% Q0 = 0.075 / 86400 * 1.1770e+04 * 40e3;
% t = 9 * 86400;
% Eprime = 6.666667e9;
% nu=0.25;
% etaprime = 2000;
% H = 5000;
% x=linspace(0,1,100);
% [L, w, p, dpdx, v_x] = PKNConstantFlux(Q0, t, Eprime, nu, etaprime, H, x);
% x=x.*L;

% Compute reduced parameters
E = Eprime * (1 - nu^2);
Ebar = (2/pi) * E / (1 - nu^2);
mubar = pi^2 * (etaprime / 12);

% Compute similarity scaling
W = ((mubar / Ebar) * (Q0^2) / H)^(1/5) * t^(1/5);  % Width scale
L = ((Ebar / mubar) * (Q0^3) / H^4)^(1/5) * t^(4/5); % Half-length scale
v_x = (4 * ((Ebar * Q0^3) / (H^4 * mubar))^(1/5)) / (5 * t^(1/5)); % Velocity

% Spatial grid
%xi = linspace(0, 1, N); % Dimensionless coordinate
xi=x;

% Compute dimensionless width profile
wBar_m = ((12/5)^(1/3)) .* (1 - xi).^(1/3) .* ...
    (1 - (1 - xi)/96 + (23*(1 - xi).^2)/64512 - (7*(1 - xi).^3)/1327104);

% Compute width profile
w = W * wBar_m;

% Compute net pressure
p = (Ebar / H) * w;

% Compute pressure gradient
wBardx = gradient(wBar_m, xi * L);
dpdx = (Ebar / H) * W * wBardx;

% %Show the volumes match
% V = trapz(x, w)*H
% V2=Q0*t

end
