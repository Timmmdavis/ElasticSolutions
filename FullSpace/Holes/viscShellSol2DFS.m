
function [ur,ut,Srr,Stt,Szz,Err,Ett,Ezz]=viscShellSol2DFS(R1,R2,eta,mu,t,rIn,p,nu)
% viscShellSol2DFS Calculate viscoelastic shell deformation and stresses
% around a cylindrical magma chamber in plane strain condition
%
% This function implements the analytical solution for a pressurized cylindrical
% magma chamber surrounded by a viscoelastic shell under plane strain conditions.
% Derived by Tim Davis (University of Bristol, 2024) following the methodology
% of Dragoni and Magnanensi (1989) but extended to 2D plane strain.
%
% The solution assumes:
% - Cylindrical geometry (infinite in z-direction)
% - Maxwell viscoelastic shell
% - Instantaneous pressurization
% - Full space solution
% - Plane strain conditions (εzz = 0)
%
% Arguments: (input)
%   R1      - Inner radius (magma chamber radius) [m]
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
%   ut      - Tangential displacement [m]
%             (Zero due to axial symmetry)
%   Srr     - Radial stress [Pa]
%             (-p inside magma chamber)
%   Stt     - Tangential stress [Pa]
%             (-p inside magma chamber)
%   Szz     - Out-of-plane stress [Pa]
%             (ν(Srr + Stt) for plane strain)
%   Err     - Radial strain [-]
%   Ett     - Tangential strain [-]
%   Ezz     - Out-of-plane strain [-]
%             (Zero for plane strain)
%
% Regions:
%   r < R1  : Magma chamber (fluid)
%   R1<r<R2 : Viscoelastic shell
%   r > R2  : Elastic host rock
%
% Notes:
%   - Solution handles both instantaneous elastic and time-dependent viscous responses
%   - Plane strain assumption means εzz = 0 and σzz = ν(σrr + σtt)
%   - Values inside magma chamber (r < R1) return NaN for displacements/strains
%   - Solution derived and validated for geological applications
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
%   [ur,ut,Srr,Stt,Szz,Err,Ett,Ezz] = viscShellSol2DFS(R1,R2,eta,mu,t,r,p,nu);
%
% References:
%   Dragoni, M., & Magnanensi, C. (1989). Displacement and stress produced by
%   a pressurized, spherical magma chamber, surrounded by a viscoelastic
%   shell. Physics of the Earth and Planetary Interiors, 56(3-4), 316-328.
%
%   Derived by Davis, T. (2024). University of Bristol.
%
% See also: viscShellSol3DFS


% Calculate bulk modulus from shear modulus and Poisson's ratio
K=((2*mu)*(1+nu))/(3*(1-(2*nu)));
InShell=rIn>R1 & rIn<R2;
ur=zeros(size(rIn));
ut=zeros(size(rIn));
Srr=zeros(size(rIn));
Stt=zeros(size(rIn));
Szz=zeros(size(rIn));
Err=zeros(size(rIn));
Ett=zeros(size(rIn));
Ezz=zeros(size(rIn));
for i=1:numel(rIn)
    r=rIn(i);
        
    q1=(2*eta);
    p1=(2*eta)/(2*mu);
    

    if InShell(i)==1  
       
        % Solution inside viscoelastic shell (R1 < r < R2)
        ur(i) =(p*(K*R2^2 + R2^2*mu - mu*r^2))/(2*K*mu*r) + (p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*((6*K^2*R1^2*R2^2 - 6*K^2*R2^4 + 6*K*R1^2*R2^2*mu - 12*K*R2^4*mu + 6*K*R2^2*mu*r^2 - 2*R1^2*R2^2*mu^2 + 2*R1^2*mu^2*r^2 - 6*R2^4*mu^2 + 6*R2^2*mu^2*r^2)/(6*p1*K^2*R1^2*R2^2 - 6*p1*K^2*R2^4 + 6*p1*K*R1^2*R2^2*mu + q1*K*R1^2*R2^2 - 12*p1*K*R2^4*mu - q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 - 2*p1*R1^2*R2^2*mu^2 + q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 - 6*p1*R2^4*mu^2 - q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2) - (3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1)/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))/(q1*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(6*p1*K^2*R1^2*R2^2 - 6*p1*K^2*R2^4 + 6*p1*K*R1^2*R2^2*mu + q1*K*R1^2*R2^2 - 12*p1*K*R2^4*mu - q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 - 2*p1*R1^2*R2^2*mu^2 + q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 - 6*p1*R2^4*mu^2 - q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2))/(2*K*mu*r*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
        Srr(i) =- p - (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - r^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
        Stt(i) =((R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(R1^2 + r^2)*(q1 + 6*K*p1 + 6*mu*p1))/(r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2)) - p);
        
        Err(i) =(p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*((- 6*K^2*R1^2*R2^2 + 6*K^2*R2^4 - 6*K*R1^2*R2^2*mu + 12*K*R2^4*mu + 6*K*R2^2*mu*r^2 + 2*R1^2*R2^2*mu^2 + 2*R1^2*mu^2*r^2 + 6*R2^4*mu^2 + 6*R2^2*mu^2*r^2)/(- 6*p1*K^2*R1^2*R2^2 + 6*p1*K^2*R2^4 - 6*p1*K*R1^2*R2^2*mu - q1*K*R1^2*R2^2 + 12*p1*K*R2^4*mu + q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 + 2*p1*R1^2*R2^2*mu^2 - q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 + 6*p1*R2^4*mu^2 + q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2) - (3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1)/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))/(q1*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(- 6*p1*K^2*R1^2*R2^2 + 6*p1*K^2*R2^4 - 6*p1*K*R1^2*R2^2*mu - q1*K*R1^2*R2^2 + 12*p1*K*R2^4*mu + q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 + 2*p1*R1^2*R2^2*mu^2 - q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 + 6*p1*R2^4*mu^2 + q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2))/(2*K*mu*r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2)) - (p*(K*R2^2 + R2^2*mu + mu*r^2))/(2*K*mu*r^2);
        Ett(i) =(p*(K*R2^2 + R2^2*mu - mu*r^2))/(2*K*mu*r^2) + (p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*((6*K^2*R1^2*R2^2 - 6*K^2*R2^4 + 6*K*R1^2*R2^2*mu - 12*K*R2^4*mu + 6*K*R2^2*mu*r^2 - 2*R1^2*R2^2*mu^2 + 2*R1^2*mu^2*r^2 - 6*R2^4*mu^2 + 6*R2^2*mu^2*r^2)/(6*p1*K^2*R1^2*R2^2 - 6*p1*K^2*R2^4 + 6*p1*K*R1^2*R2^2*mu + q1*K*R1^2*R2^2 - 12*p1*K*R2^4*mu - q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 - 2*p1*R1^2*R2^2*mu^2 + q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 - 6*p1*R2^4*mu^2 - q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2) - (3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1)/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))/(q1*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(6*p1*K^2*R1^2*R2^2 - 6*p1*K^2*R2^4 + 6*p1*K*R1^2*R2^2*mu + q1*K*R1^2*R2^2 - 12*p1*K*R2^4*mu - q1*K*R2^4 + 6*p1*K*R2^2*mu*r^2 - 2*p1*R1^2*R2^2*mu^2 + q1*R1^2*R2^2*mu + 2*p1*R1^2*mu^2*r^2 - q1*R1^2*mu*r^2 - 6*p1*R2^4*mu^2 - q1*R2^4*mu + 6*p1*R2^2*mu^2*r^2 + q1*R2^2*mu*r^2))/(2*K*mu*r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));


    else
        % Solution in elastic host rock (r > R2)
        ur(i) =(R2^2*p)/(2*mu*r) + (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - R2^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(2*mu*r*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
        Srr(i) =- (R2^2*p)/r^2 - (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - R2^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
        Stt(i) =(R2^2*p)/r^2 + (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - R2^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));

        Err(i) =- (R2^2*p)/(2*mu*r^2) - (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - R2^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(2*mu*r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
        Ett(i) =(R2^2*p)/(2*mu*r^2) + (R2^2*p*q1*exp(-(t*(3*K*R2^2*q1 - 3*K*R1^2*q1 + R1^2*mu*q1 + 3*R2^2*mu*q1 + 12*K*R1^2*mu*p1))/(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2))*(R1^2 - R2^2)*(cosh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1))) + (sinh((q1*t*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2))/(2*mu*p1*(R1^2*q1 + 3*R2^2*q1 + 6*K*R1^2*p1) - R1^2*q1*(q1 + 6*K*p1) + R2^2*q1*(q1 + 6*K*p1)))*(- 18*p1*K^2*R1^2 + 18*p1*K^2*R2^2 - 24*p1*K*R1^2*mu - 3*q1*K*R1^2 + 36*p1*K*R2^2*mu + 3*q1*K*R2^2 + 6*p1*R1^2*mu^2 - 7*q1*R1^2*mu + 18*p1*R2^2*mu^2 + 3*q1*R2^2*mu))/((q1 + 6*K*p1 + 6*mu*p1)*(9*K^2*R1^4 - 18*K^2*R1^2*R2^2 + 9*K^2*R2^4 + 6*K*R1^4*mu - 24*K*R1^2*R2^2*mu + 18*K*R2^4*mu + R1^4*mu^2 + 6*R1^2*R2^2*mu^2 + 9*R2^4*mu^2)^(1/2)))*(q1 + 6*K*p1 + 6*mu*p1))/(2*mu*r^2*(R2^2*q1^2 - R1^2*q1^2 - 6*K*R1^2*p1*q1 + 6*K*R2^2*p1*q1 + 2*R1^2*mu*p1*q1 + 6*R2^2*mu*p1*q1 + 12*K*R1^2*mu*p1^2));
       
    end

    Szz(i)=nu*(Srr(i)+Stt(i));%Pollard 8.35
    Ezz(i)=0; %Plane strain    

end
%In the fluid
ur(rIn<R1)=nan;
Srr(rIn<R1)=-p;
Stt(rIn<R1)=-p;
Err(rIn<R1) =nan;
Ett(rIn<R1) =nan;

end
