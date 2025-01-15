function [uz, ur] = PennyNearSurfacePlate(d, a, nu, E, p, r)
% PennyNearSurfacePlate Calculate displacements due to a near-surface 
% penny-shaped hydraulic fracture using the plate theory approximation
%
% This function implements the asymptotic solution for a near-surface 
% penny-shaped hydraulic fracture as described in:
% "An Asymptotic Solution for the Problem of a Hydraulic Fracture in a 
% Elastic Layer" by Andrew P. Bunger and Emmanuel Detournay
%
% The solution approximates the elastic layer as a plate, valid when the 
% fracture radius is much larger than the layer thickness (a/d > 5)
%
% Arguments: (input)
%   d       - Depth of fracture center from surface [m]
%   a       - Radius of penny-shaped fracture [m]
%   nu      - Poisson's ratio of medium [-]
%   E       - Young's modulus [Pa]
%   p       - Change in hydrostatic pressure within fracture [Pa]
%   r       - Radial distance from center where displacement is calculated [m]
%
% Arguments: (output)
%   uz      - Vertical displacement at surface [m]
%   ur      - Radial displacement at surface [m]
%
% Notes:
%   - Solution is based on plate theory approximation
%   - Valid when a/d > 5 (fracture radius much larger than depth)
%   - Assumes uniform pressure distribution in fracture
%   - Uses small strain theory
%
% Example usage:
%   d = 100;             % Depth of 100m
%   a = 1000;           % Radius of 1000m
%   nu = 0.25;          % Poisson's ratio
%   E = 30e9;           % Young's modulus of 30 GPa
%   p = 1e6;            % Pressure change of 1 MPa
%   r = linspace(0,2000,100); % Observation points
%   [uz,ur] = PennyNearSurfacePlate(d, a, nu, E, p, r);
%   plot(r, uz);        % Plot vertical displacement profile
%
% References:
%   Bunger, A.P., and Detournay, E., "An Asymptotic Solution for the 
%   Problem of a Hydraulic Fracture in a Elastic Layer"
%
% See also: None

    % Rename variables to match paper notation
    H = d;      % Layer thickness/depth
    R = a;      % Fracture radius
    
    % Calculate non-dimensional radial coordinate (Eq. not numbered)
    rho = r./R;
    
    % Calculate plate flexural rigidity (Eq. 1)
    D = (E/(1-nu^2)*H^3)/12;
    
    % Empirical coefficient from numerical calculations (Eq. 15)
    delta = 0.620;
    
    % Calculate dimensionless parameter η (Eq. 21)
    eta = (4.*delta)./(3+2.*nu).*H./R;
    
    % Calculate vertical displacement (Eq. 22)
    uz = (R.^4.*p)./(64.*D).*((1-rho.^2).^2 + 4.*eta.*(1-rho.^2));
    
    % Calculate radial displacement from vertical displacement gradient
    % First calculate dw/dr
    dwdr = (R^4.*p.*((4.*r.*(r.^2./R.^2 - 1))./R.^2 - ...
            (8.*eta.*r)./R.^2))./(64.*D);
    
    % Calculate radial displacement based on plate theory
    ur = -(dwdr/2*d);
end