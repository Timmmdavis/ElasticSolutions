function [ur, Srr, Stt] = viscShellSol3DFS(R1, R2, eta, mu, t, rIn, p, nu)
% viscShellSol3DFS Calculate viscoelastic shell deformation around a pressurized 
% spherical magma chamber in a full space
%
% This function implements the analytical solution for a pressurized spherical
% magma chamber surrounded by a viscoelastic shell, as described in:
% "Earthquake and Volcano Deformation" (Segall, 2010) and
% "Displacement and stress produced by a pressurized, spherical magma chamber, 
% surrounded by a viscoelastic shell" (Dragoni and Magnanensi, 1989)
%
% The solution assumes:
% - Spherically symmetric geometry
% - Maxwell viscoelastic shell
% - Instantaneous pressurization
% - Full space solution
%
% Arguments: (input)
%   R1      - Radius of magma chamber [m]
%   R2      - Outer radius of viscoelastic shell [m]
%   eta     - Viscosity of shell [Pa·s]
%   mu      - Shear modulus of shell [Pa]
%   t       - Time elapsed since pressurization [s]
%   rIn     - Radial distances for calculation points [m]
%   p       - Magnitude of pressure change [Pa]
%   nu      - Poisson's ratio of shell [-]
%
% Arguments: (output)
%   ur      - Radial displacement [m]
%             (NaN inside magma chamber, r < R1)
%   Srr     - Radial stress [Pa]
%             (-p inside magma chamber)
%   Stt     - Tangential stress (σφφ = σθθ) [Pa]
%             (-p inside magma chamber)
%
% Regions:
%   r < R1  : Magma chamber (fluid)
%   R1<r<R2 : Viscoelastic shell
%   r > R2  : Elastic host rock
%
% Notes:
%   - Solution includes both instantaneous elastic and time-dependent viscous responses
%   - Tangential displacement is zero due to spherical symmetry
%   - K (bulk modulus) is calculated from mu and nu
%   - Values inside magma chamber (r < R1) return NaN for displacement
%
% Example usage:
%   R1 = 1000;           % 1 km chamber radius
%   R2 = 2000;           % 2 km shell radius
%   eta = 1e18;          % 10^18 Pa·s viscosity
%   mu = 10e9;           % 10 GPa shear modulus
%   t = 365*24*3600;     % 1 year in seconds
%   nu = 0.25;           % Poisson's ratio
%   p = 10e6;            % 10 MPa pressure change
%   r = linspace(0,5000,100); % Observation points
%   [ur,Srr,Stt] = viscShellSol3DFS(R1,R2,eta,mu,t,r,p,nu);
%
% References:
%   Segall, P. (2010). Earthquake and Volcano Deformation. 
%   Princeton University Press.
%
%   Dragoni, M., & Magnanensi, C. (1989). Displacement and stress produced by
%   a pressurized, spherical magma chamber, surrounded by a viscoelastic
%   shell. Physics of the Earth and Planetary Interiors, 56(3-4), 316-328.
%
% See also: None

    % Calculate bulk modulus from shear modulus and Poisson's ratio
    K1 = ((2*mu)*(1+nu))/(3*(1-(2*nu)));
    K = K1;
    
    % Initialize output arrays
    InShell = rIn>R1 & rIn<R2;
    ur = zeros(size(rIn));
    Srr = zeros(size(rIn));
    Stt = zeros(size(rIn));
    
    % Calculate solutions for each radial position
    for i = 1:numel(rIn)
        r = rIn(i);
        if InShell(i) == 1  
            % Solution inside viscoelastic shell (R1 < r < R2)
            ur(i) = (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))* ...
                    (4*mu*p*r^3 + 3*K*R1^3*p - 3*K*R2^3*p - 4*R2^3*mu*p))/ ...
                    (12*K*mu*r^2) + (p*(3*K*R2^3 + 4*R2^3*mu - 4*mu*r^3))/(12*K*mu*r^2);
            
            Srr(i) = -p - (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))* ...
                     (p*R1^3 - p*r^3))/r^3;
            
            Stt(i) = (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))* ...
                     (p*R1^3 + 2*p*r^3))/(2*r^3) - p;
        else
            % Solution in elastic host rock (r > R2)
            ur(i) = (R2^3*p)/(4*mu*r^2) + (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + ...
                    4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/(4*mu*r^2);
            
            Srr(i) = -(R2^3*p)/r^3 - (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + ...
                     4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/r^3;
            
            Stt(i) = (R2^3*p)/(2*r^3) + (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + ...
                     4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/(2*r^3);
        end
    end
    
    % Set values inside magma chamber (r < R1)
    ur(rIn<R1) = nan;
    Srr(rIn<R1) = -p;
    Stt(rIn<R1) = -p;
end