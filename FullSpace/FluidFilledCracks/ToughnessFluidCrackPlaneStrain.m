function [a, wf, dwfdt, dpdx, vx] = ToughnessFluidCrackPlaneStrain(q0, t, Kc, Eprime, etaprime, x)
% ToughnessFluidCrackPlaneStrain: Computes crack half-length, fluid width, width rate of change, 
% and pressure gradient in a 2D plane strain fracture under toughness-dominated conditions driven 
% by a constant flux.
%
%    Diagram: cross-section through the crack
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
%     Kc        - Fracture toughness (Pa·sqrt(m)).
%     Eprime    - Plane strain Young's modulus (E/(1-nu^2)) (Pa).
%     etaprime  - Viscosity parameter for pressure gradient computation (12*eta) (Pa·s).
%     x         - Spatial coordinate along the crack (-1 to 1).
% 
% Arguments: (output)
%     a         - Half-length of the crack (m).
%     wf        - Fluid width (separation between walls) at given time (m).
%     dwfdt     - Rate of change of width with time.
%     dpdx      - Pressure gradient in the crack (Pa/m).
%     vx        - Tip speed (m/s).
%
% References:
% Tada, (2000): Stress analysis of cracks handbook.
% 
% Example usage:
% q0 = 1e-6;  % m^2/s
% t = 10;     % s
% Kc = 1e6;   % Toughness
% Eprime = 1e9; % Plane strain Young's modulus
% etaprime = 0.001;
% x = linspace(-1,1,100); % Position array
% [a, wf, dwfdt, dpdx, vx] = PlaneStrainCrack(q0, t, Kc, Eprime, etaprime, x);
% x=x*a;
%
% Author: Tim Davis

% Compute half-length of the crack
a = ((Eprime .* q0 .* t) / (Kc .* pi.^(1/2))).^(2/3);
x=x*a;

%Tip speed - diff(a)./diff(t)
vx=(2*Eprime*q0)/(3*Kc*pis^(1/2)*((Eprime*q0*t)/(Kc*pis^(1/2)))^(1/3));

% Compute fluid width
wf = (4 .* (Kc).^(4/3)) ./ ((pi .* q0 .* t).^(1/3) .* Eprime.^(4/3)) .* sqrt(a.^2 - x.^2);

% Compute rate of change of width with time
dwfdt = 4 .* q0 .* (a.^2 + x.^2) ./ (3 .* pi .* a^2 .* sqrt(a.^2 - x.^2));

% Compute pressure gradient
dpdx = etaprime .* dwfdt ./ wf.^3;

end
