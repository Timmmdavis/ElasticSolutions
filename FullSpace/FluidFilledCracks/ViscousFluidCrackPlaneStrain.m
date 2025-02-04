function [a, wf, p, q, dpdx, vx] = ViscousFluidCrackPlaneStrain(q0, t, Eprime, etaprime, x)
% ViscousDominatedCrack: Computes crack half-length, fluid width, pressure, 
% flux distribution, and pressure gradient in a 2D plane strain fracture 
% dominated by viscous forces driven by a constant flux. Note this is only the
% first-order solution. 
% 
%     Diagram: cross-section through the crack
%
%                         q₀ (fluid flux)
%                            ↓
%                    Injection point (x = 0)
%                            *
%                            │
%                     ────────────────   
%                ────        │         ────                    
%      vx ←  ────         wf │             ────  ← top fracture wall
%         ←  ────            │             ────  ← bottom fracture wall      
%               ────         │         ────           
%                    ─────────────────
%                            *
%                            x ───────────────→
%                                     a
%                              (crack half–length)
% 
%          Pressure gradient: dp/dx along the crack (x-axis)
%          Tip speed at fracture tips: vx (at x = ±a)
%
% Arguments: (input)
%     q0        - Constant fluid flux into one wing (m^2/s).
%     t         - Time (s).
%     Eprime    - Plane strain Young's modulus (E/(1-nu^2)).
%     etaprime  - Viscosity parameter for pressure gradient computation (eta*12).
%     x         - Spatial coordinate along the crack (0 to l).
% 
% Arguments: (output)
%     a         - Half-length of the crack (m).
%     wf        - Fluid width (separation between walls) at given time (m).
%     p         - Pressure distribution along the crack (Pa).
%     q         - Flux distribution along the crack (m^2/s).
%     dpdx      - Pressure gradient in the crack (Pa/m)..
%     vx        - Tip speed (m/s).
% 
% References:
% Adachi, J. I. (2001) - Fluid-driven fracture in permeable rock. https://bris.idm.oclc.org/login?url=https://www.proquest.com/dissertations-theses/fluid-driven-fracture-permeable-rock/docview/252091478/se-2 
% Detournay (2004) - https://doi.org/10.1061/(ASCE)1532-3641(2004)4:1(35)
% Zolfaghari, Dontsov and Bunger (2017) -  https://doi.org/10.1002/nag.2755
%
% Example usage:
% q0 = 1e-6;  % m^2/s
% t = 10;     % s
% Eprime = 1e9; % Plane strain Young's modulus
% etaprime = 0.001;
% x = linspace(0,1,100); % Position array
% [l, w, p, q, dpdx, vx] = ViscousFluidCrackPlaneStrain(q0, t, Eprime, etaprime, x);
% x=x*a;
%
% Author: Tim Davis

% Constants
M = -0.15601;
N = 0.066322;

% Calculate l(t)
a =  0.6162 .* ((Eprime .* q0.^3) ./ etaprime).^(1./6) .* t.^(2./3);
x=x*a;

%Tip speed diff(a)/diff(t)
vx=(2*0.6162 *((Eprime*q0^3)/etaprime)^(1/6))/(3*t^(1/3));

lt = a;
term1 = 0.6162 .* ((etaprime .* q0.^3) ./ Eprime).^(1./6) .* t.^(1./3);
term2 = sqrt(3) .* (1 - (x./lt).^2).^(2./3);
term3 = M .* (1 - (x./lt).^2).^(5./3);
term4 = N .* (4 .* sqrt(1 - (x./lt).^2) + 2.*(x./lt).^2 .* ...
        log((1 - sqrt(1 - (x./lt).^2))./(1 + sqrt(1 - (x./lt).^2))));
wf = term1 .* (term2 + term3 + term4);

term1 = ((etaprime .* Eprime.^2 )./ t).^(1./3);
term2 = 1./(3.*pi) .* beta(1./2, 2./3);
term3 = sqrt(3) .* hypergeom([-1./6, 1], 1./2, (x./lt).^2);
term4 = (10./7) .* M .* hypergeom([-7./6, 1], 1./2, (x./lt).^2);
term5 = N .* (2 - pi .* abs(x./lt));
p = term1 .* (term2 .* (term3 + term4) + term5);

term1 = 0.3797 .* q0;
term2 = (2./7) .* (sqrt(3) + 10.*M./13) .* beta(1./2, 2./3);
term3 = -sqrt(3) .* (x./lt) .* hypergeom([1./2, -2./3], 3./2, (x./lt).^2);
term4 = (2./sqrt(3)) .* (x./lt) .* (1 - (x./lt).^2).^(2./3);
term5 = -M .* (x./lt) .* hypergeom([1./2, -5./3], 3./2, (x./lt).^2);
term6 = (2.*M./3) .* (x./lt) .* (1 - (x./lt).^2).^(5./3);
term7 = (4.*N./3) .* (acos(x./lt) + (x.^3)./(2.*lt.^3) .* ...
        log((1 - sqrt(1 - (x./lt).^2))./(1 + sqrt(1 - (x./lt).^2))));
q = term1 .* (term2 + term3 + term4 + term5 + term6 + term7);

dpdx=gradient(p)./gradient(x);

end