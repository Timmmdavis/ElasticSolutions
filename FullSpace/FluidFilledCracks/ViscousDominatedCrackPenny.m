 function [tmk, V, a, vr, K, P, wf, dpdr] = ViscousDominatedCrackPenny(Q, t, Eprime, etaprime, E, nu, r, Kprime)
% ViscousDominatedCrackPenny: Computes fracture metrics for a viscous-dominated,
% penny-shaped fracture driven by a constant flux.
%
% This function determines the transition time scale for a toughness-dominated
% fracture regime, calculates the fracture volume using both pressure-based and
% width-based integration methods, and computes key fracture parameters in the
% viscous-dominated regime following Detournay (2016). In particular, it computes
% the fracture radius, propagation velocity, stress intensity factor (K), and the
% pressure and width (aperture) distributions.
%
%
%   Diagram: cross-section through the penny shaped crack
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
%
% Inputs:
%   Q       - Injection rate (m^3/s)
%   t       - Time (s)
%   Eprime  - Plane strain Young's modulus (E/(1-nu^2)) (Pa)
%   etaprime -Fluid viscosity (12*eta) (Pa·s)
%   E       - Young's modulus (Pa)
%   nu      - Poisson's ratio (dimensionless)#
%   r       - Radial coordinate within the fracture (0 to 1).
%   Kprime  - Fracture toughness parameter (Kprime=4*sqrt(2/pi)*Kc) (Pa·sqrt(m))
%
% Outputs:
%   tmk     - Toughness-dominated transition time scale (s)
%   V       - Injected fracture volume calculated as Q*t (m^3)
%   R       - Fracture radius (m)
%   vr      - Fracture (radial) propagation velocity (m/s)
%   K       - Stress intensity factor computed via integration (Pa·sqrt(m))
%   P       - Pressure distribution within the fracture (Pa)
%   wf      - Fracture aperture (fluid width) distribution (m)
%
% References:
% Detournay (2016) - https://doi.org/10.1146/annurev-fluid-010814-014736
% Savitski and Detournay (2001) - https://doi.org/10.1016/S1620-7742(01)01323-X
% Savitski and Detournay (2001) - https://doi.org/10.1016/S0020-7683(02)00492-4
%
% Example usage:
% Q       = 1e-3;            % m^3/s
% t       = 100;             % s
% Kprime  = 1e6;             % Pa·m^(0.5)
% E       = 1e9;             % Pa
% nu      = 0.25;
% Eprime  = E/(1 - nu^2);      % Pa
% etaprime = 0.001;           % Pa·s
% r=linspace(0,1,100)
% [tmk, V, a, vr, K, P, wf] = ViscousDominatedCrackPenny(Q, t, Eprime, etaprime, E, nu, r, Kprime);
% r=r*a;
%
% Author: Tim Davis
   
if nargin<7
    Kprime=0;
end
Kc=Kprime/(4*sqrt(2/pi));

tmk=sqrt((Eprime^13*etaprime^5*Q^3)/(Kprime^(18)));

if t>tmk
    disp('Use toughness dominated regime')
end    

 %% Volume calculations
%Internal pressure from given area
%Davis 2022 JFM manuscript, penny-crack
V=Q*t;  %Volume 

%   VNum    - Fracture volume computed via integration of the pressure-based formulation (m^3)
%   VNum2   - Fracture volume computed via integration of the width-based formulation (m^3)

% %Using the ring dipole in the tada book show the pressure is correct
% %such that we can calcuate the correct volume.
% fun = @(r) CalculateVolume(r,R,etaprime,Eprime,t,Q,mscl,E,nu);
% VNum=quadgk(fun,0,R); 
% %Using the widths from the Savitski 2001 paper to calc volume - shell
% %integration
% fun = @(r) r.*CalculateWidth(r,R,etaprime,Eprime,t,Q,mscl);
% VNum2=(2*pi)*quadgk(fun,0,R);    
% %Quick test to make sure the volumes are all the same. 
% if isnan(VNum)
%     if round(V/VNum)~=1 && round(V/VNum2)~=1
%         error('Both pressure and widths are incorrect now')
%     elseif round(V/VNum)~=1 
%         error('Pressure is incorrect now')        
%     elseif round(V/VNum2)~=1
%         error('Width is incorrect now')
%     end
% end

%% Viscous dominated Detournay 2016
mscl=0.6955; %Eq.35
a=mscl*((Eprime*Q^3*t^4)/(etaprime))^(1/9);%Eq.31 - fracture radius
[vr]=CalculateFrontVelocity(mscl,Eprime,Q,etaprime,t);
r=r*a;

%% K calculations
%Using the ring dipole in the Tada book show the analytical pressure
%gives correct K: should be approx 0.
fun = @(r) CalculateK(r,a,etaprime,Eprime,t,Q,mscl,E,nu);
K=quadgk(fun,0,a); 

%Using the widths from the Savitski 2001 paper to calc volume - shell
%integration
[P,wf] = CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);

%Calculate the pressure gradient across the fracture.
[dpdr]=CalculateDpdr(r,a,etaprime,Eprime,t,Q,mscl);

% %% Some Energy analysis
% % Use pressure as close to src as possible instead!!!
% rho=eps();
% w1=2.479;   A1=3.581*10^-1; B1=9.269*10^-2;
% muPscl=A1*(w1-(2/(3*(1-rho)^(1/3))))-B1*(log(rho/2)+1);
% P3=((Eprime^2*etaprime)/t)^(1/3)*muPscl; %From Mori table 3 -3-D buoyant hydraulic fractures: constant release
% Wp=Q*abs(P3);
% % Elastic work - using Sneddon integration technique - you need two time steps.
% fun = @(r) r.*CalculateElasticEnergy(r,a,etaprime,Eprime,t,Q,mscl,E,nu,vr);
% Wet1=(2*pi*quadgk(fun,0,a));
% dt=t/100; %go some finite time forwards...
% fun = @(r) r.*CalculateElasticEnergy(r,a,etaprime,Eprime,t+dt,Q,mscl,E,nu,CalculateFrontVelocity(mscl,Eprime,Q,etaprime,t+dt));
% Wet2=(2*pi*quadgk(fun,0,a));
% We=(Wet2-Wet1)/dt;
% %Fluid work
% fun = @(r) r.*CalculateDissipation(r,a,etaprime,Eprime,t,Q,mscl);
% Wmu=2*pi*quadgk(fun,0,a);
% %Fracture energy
% Wk=(2*pi*a)*(Kc^2/Eprime)*vr;


end

function [vr]=CalculateFrontVelocity(mscl,Eprime,Q,etaprime,t)
    vr=mscl*4/9*((Eprime*Q^3)/(etaprime*t^5))^(1/9);%diff(R)/diff(t)
end

function [P,wf]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl)


    %% Savitski 2002 "if not stated where from" && Savitski 2001 - "if it states: C. R. Acad. Sci. Paris,"
    %Table 1
    w1=2.479; 
    %Constants - Table 3    
    A1=3.581*10^-1;
    B1=9.269*10^-2;
    C1=6.846*10^-1;
    C2=7.098*10^-2;
    %table.4 (diff to other one in det 2016)
    %mscl=6.955*10^-1; 

    %Below eq.14
    rhoV=r./a;
    %Begin loop
    muPscl=zeros(size(rhoV));
    muWscl=zeros(size(rhoV));
    for j=1:numel(rhoV) 
        rho=rhoV(j);
        if rho>1-0.1%Detournay 2006 - "noting that the m-asymptote applies to approximately 10% of the fracture radius from the tip" 
         %If you want to use tip sol
            %Asymptotic sols of DETOURNAY 2016
            muPscl(j)=-3^(-4/3)*(1-rho)^(-1/3);
            muWscl(j)=2*3^(1/6)*mscl*(1-rho)^(2/3); %Eq.34 Detournay 2016 - approx close 2 tip
        else 
            %Eq.69 / (Eq.30 in C. R. Acad. Sci. Paris, t. 329, Série II b,
            %p. 255–262, 2001)
            muPscl(j)=A1*(w1-(2/(3*(1-rho)^(1/3))))-B1*(log(rho/2)+1);
            %Eq.68 
            muWscl(j)=(sqrt(70)/3*C1+(4*sqrt(5))/9*C2*(13*rho-6))*(1-rho)^(2/3)+B1*(8/pi*sqrt(1-rho)-8/pi*rho*acos(rho));
            %Eq.31 - C. R. Acad. Sci. Paris, t. 329, Série II b, p. 255–262, 2001
            muWscl(j)=(sqrt(70)/3*C1+(4*sqrt(5))/9*C2*(13*rho-6))*(1-rho)^(2/3)+B1*(-4*rho+8/pi*sqrt(1-rho^2)+8/pi*rho*atan(rho/sqrt(1-rho^2)));
        end
        if rho>1-0.9 %Singular pressure at source - See Savitski 2002 eq.48
            muPscl(j)=-log(rho);
            % C11=6.846*1e-1;
            % C21=7.098*1e-2;
            % B11=9.269*1e-2;
            % muWscl(j)=(sqrt(70)/3*C11+(4*sqrt(5))/9*C21*(13*rho-6))*(1-rho)^(2/3)+B11*((8/pi)*(1-rho)^(1/2)-(8/pi)*rho*acos(rho));
        end

    end    

%     %Eq.27 
%     em=(etaprime/(Eprime*t))^(1/3);
%     Lm=((Eprime*Q^3*t^4)/(etaprime))^(1/9);
%
%     %Eq.12:14
%     P=em.*Eprime.*muPscl;
%     wf=em.*Lm.*muWscl;  


%     T=etaprime/Eprime;%Eq.7 - C. R. Acad. Sci. Paris, t. 329
%     tau=t/T;
    tau=t*Eprime/etaprime;
    %Eq.10 - C. R. Acad. Sci. Paris, t. 329, Série II b, p. 255–262, 2001
    %P=Eprime.*tau^-(1/3)*muPscl;
    P=(Eprime.*(A1.*(w1 - 2./(3*(1 - r./a).^(1/3))) - B1.*(log(r./(2.*a)) + 1)))./((Eprime.*t)./etaprime).^(1/3);

%     gamma=6.955.*10^-1;%Savitski 2002 - table 4
%     L=R./(gamma.*tau.^(1/3+1/9));%eq 11
%     wf=L.*gamma.*tau.^(1/3-2/9).*muWscl;   
    %Simplified:
    wf=a*tau^(-1/3)*muWscl; 
   
end

function [dpdr]=CalculateDpdr(r,a,etaprime,Eprime,t,Q,mscl)

    [~,wf]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);

    %Analytical diff: to get grad P. You could use: "dpdr=gradient(P)./gradient(r);"
    %Constants - Table 3    
    A1=3.581*10^-1;
    B1=9.269*10^-2;
    dpdr=-(Eprime.*(B1./r + (2.*A1)./(9.*a.*(1 - r./a).^(4/3))))./((Eprime.*t)./etaprime).^(1/3);

end

function [dissipation]=CalculateDissipation(r,a,etaprime,Eprime,t,Q,mscl)

    [~,wf]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);
    [dpdr]=CalculateDpdr(r,a,etaprime,Eprime,t,Q,mscl);

    %Rich's way
    dissipation=(dpdr.^2.*wf.^3)./(etaprime);

    %My way...
    q=(dpdr.*wf.^3)./(etaprime);
    % q2=Q./(2*pi*r); %Close to well bore - A.P. Bunger / International Journal of Solids and Structures 50 (2013) 1538–1549 eq.A14 
    % q(r<R/15)=q2(r<R/15); 
    dissipation=(q.^2.*etaprime)./wf.^3;

end

function [K]=CalculateK(r,a,etaprime,Eprime,t,Q,mscl,E,nu)
    
    %Pressure from simple first order analytical solution
    [P,~]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);
    %Using this pressure to calculate K via integration
    %Ring dipole: Tada 2000 p.346 24.5.
    K=(2.*P)./sqrt(pi.*a).*((r)./(sqrt(a.^2-r.^2)));
   
end

function [V]=CalculateVolume(r,a,etaprime,Eprime,t,Q,mscl,E,nu)
    
    %Pressure from simple first order analytical solution
    [P,~]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);
    %Using this pressure to calculate the volume via integration
    %Ring dipole: Tada 2000 p.346 24.5.
    V=((16.*(1-nu.^2))./E).*P.*r.*sqrt(a.^2-r.^2); 
   
end

function [wf]=CalculateWidth(r,a,etaprime,Eprime,t,Q,mscl)
    
    [~,wf]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);
   
end

function [We]=CalculateElasticEnergy(r,a,etaprime,Eprime,t,Q,mscl,E,nu,vy)
    
    %Pressure from simple first order analytical solution
    [P,wf]=CalculatePressureAndWidth(r,a,etaprime,Eprime,t,Q,mscl);

    %Sneddon 1965 -sec.4 after eq.13 - shell integation of this is the
    %energy
    We=P.*wf;
   
end