function [a, wf, vy, P, dP, dpdx] = ToughnessDominatedCrackPenny(Q, t, Kprime, Eprime, nu, r, etaprime)
% ToughnessDominatedCrackPenny: Computes fracture radius, fluid width, 
% velocity, pressure, and pressure rate of change in a toughness-dominated 
% penny-shaped fracture driven by a constant flux.
% 
%
%    Diagram: cross-section through the penny shaped crack
%
%                         Q (fluid flux)
%                            ↓
%                    Injection point (r = 0)
%                            *
%                            │
%                     ────────────────   
%                ────        │         ────                    
%      vr ←  ────         wf │             ────  ← top fracture wall
%         ←  ────            │             ────  ← bottom fracture wall      
%               ────         │         ────           
%                    ─────────────────
%                            *
%                            r ───────────────→
%                                     a
%                              (crack radius)
% 
%          Pressure gradient: dp/dr along the crack (r-axis)
%          Tip speed at fracture tips: vr (at r = ±a)
%
% Arguments: (input)
%     Q        - Injection rate (m^3/s).
%     t        - Time (s).
%     Kprime   - Fracture toughness parameter. (Kprime=4*sqrt(2/pi)*Kc) (Pa·sqrt(m))
%     Eprime   - Plane strain Young's modulus (E/(1-nu^2)).
%     nu       - Poisson's ratio.
%     r        - Radial coordinate within the fracture (0 to 1).
%     etaprime - Fluid viscosity (12*eta) (Pa·s)
% 
% Arguments: (output)
%     a        - Fracture radius (m).
%     wf       - Fluid width (separation between walls) at given time (m).
%     vy       - Radial velocity of fracture propagation (m/s).
%     P        - Pressure inside the fracture (Pa).
%     dP       - Rate of change of pressure over time at source (Pa/s).
%     dpdx     - Gradient in pressure along fracture radius (Pa/m).
% 
% References:
% Tada, (2000): Stress analysis of cracks handbook.
% 
% Example usage:
% Q = 1e-3;  % m^3/s
% t = 100;   % s
% Kprime = 1e6;
% Eprime = 1e9;
% nu = 0.25;
% etaprime=100*12;
% r = linspace(0,1,100); % Radial position array
% [R, wf, vr, P, dP, dpdr] = ToughnessDominatedCrackPenny(Q, t, Kprime, Eprime, nu, r, etaprime);
% r=r*R;
%
% Author: Tim Davis


E=Eprime*(1-nu^2);
% Compute shear modulus
mu = E / (2 * (1 + nu));

% Compute fracture radius: Using Tada stress analysis of cracks handbook: Eq.36
kscl=(9/(pi^2*2))^(1/5);
a=kscl*((Eprime^2*Q^2*t^2)/(Kprime^2))^(1/5);

r=r*a; 

%Volume
V=Q*t;   

% Compute pressure inside the fracture -Davis 2022 JFM manuscript, penny-crack
P=((3*Eprime)/(16))*(V/a^3);%Tada 342

%Tada 2000 P.342 - both walls
wf=((8)./(pi.*Eprime)).*P.*sqrt(a.^2-r.^2);

% %diff(wf)/diff(t)
% dwfdt=(3.*Q.*(a.^2 - r.^2).^(1/2))./(2.*a.^3.*pi);

%dpdx = etaprime * dwfdt / wf^3;
dpdx=(4.*a.^6.*etaprime.*pi.^2)./(9.*Q.^2.*t.^3.*(a.^2 - r.^2));
   
% Compute radial velocity
vy=2/5*(9/(pi^2*2)*(Q^2*Eprime^2)/(t^3*Kprime^2))^(1/5);%diff(R)/diff(t)    

% Compute rate of change of pressure over time
dP=(3*Eprime^2*Q^3*mu*t^2)/(40*Kprime^2*kscl^3*(nu - 1)*((Eprime^2*Q^2*t^2)/Kprime^2)^(8/5));

% %% Some energy analysis:
% Wp=Q*P;  
% We=16*(1-nu^2)*P^2*(a^2*vy)/(3*E); %Sneddon -  a note on the problem of a... -> to work per s BOTH SIDES!
% Kc=Kprime/(4*sqrt(2/pi));
% Wk=(2*pi*a)*(Kc^2/Eprime)*vy;


end
