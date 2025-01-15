function [u, v, w, dV, K1, K2] = fialko2001disp(x0,y0,z0,p0,mu,a,nu,x,y,z)
% 3D Green's function for sill-like source (Fialko et al., 2001)
% all parameters are in SI (MKS) units
%
% INPUT
% x0,y0     coordinates of the center of the sphere 
% z0        depth of the center of the sill (positive downward and
%              defined as distance below the reference surface)
% p0        excess pressure 
% mu        shear modulus
% a         radius of the sphere
% nu        Poisson's ratio
% x,y       benchmark location
% z         depth within the crust (z=0 is the free surface)
% 
% OUTPUT
% u         horizontal (East component) deformation
% v         horizontal (North component) deformation
% w         vertical (Up component) deformation
% dV        volume change
% K1        Mode 1 stress intensity 
% K1        Mode 2 stress intensity
% *************************************************************************
% Fialko, Y, Khazan, Y and M. Simons (2001). Deformation due to a 
% pressurized horizontal circular crack in an elastic half-space, with 
% applications to volcano geodesy. Geophys. J. Int., 146, 181–190
% *************************************************************************
%==========================================================================
% USGS Software Disclaimer 
% The software and related documentation were developed by the U.S. 
% Geological Survey (USGS) for use by the USGS in fulfilling its mission. 
% The software can be used, copied, modified, and distributed without any 
% fee or cost. Use of appropriate credit is requested. 
%
% The USGS provides no warranty, expressed or implied, as to the correctness 
% of the furnished software or the suitability for any purpose. The software 
% has been tested, but as with any complex software, there could be undetected 
% errors. Users who find errors are requested to report them to the USGS. 
% The USGS has limited resources to assist non-USGS users; however, we make 
% an attempt to fix reported problems and help whenever possible. 
%==========================================================================
P_G=p0/mu;

% General parameters ******************************************************
eps = 1E-8;                                                                 % relative accuracy 
rd = a;                                                                     % avoid conflict with line 49 
h = z0/rd;                                                                  % dimensionless source depth                
% *************************************************************************

% Coordinates transformation **********************************************
x = (x-x0)/rd; y = (y-y0)/rd; z = (z-z0)/rd;                                % translate and scale locations
r = sqrt(x.^2+y.^2);                                                        % compute radial distance 
% *************************************************************************

% solve for PHI and PSI, Fialko et al. (2001), eq. (26) *******************
[csi1, w1] = gauleg(eps,10,41);                                                 
[csi2, w2] = gauleg(10,60,41);                                                  
csi = cat(2,csi1,csi2);  wcsi = cat(2,w1,w2);                               % ascissas and weights for Gauss-Legendre quadrature 

if size(csi,1)==1
    csi = csi'; 
end                                         % check that csi is a column vectors
[phi, psi, t, wt] = psi_phi(h);
%See appendix - these are non-dim stress intensities - Just above appendix B:
PHI = sin(csi*t)*(wt'.*phi);                                                % Gauss-Legendre quadrature
PSI = (sin(csi*t)./(csi*t) - cos(csi*t))*(wt'.*psi);                        % Gauss-Legendre quadrature
K1=-(wcsi*PHI)*p*sqrt(a*pi); %Redim...
K2=-(wcsi*PSI);
% *************************************************************************

% compute A and B, Fialko et al. (2001), eq. (24) *************************
% NOTE there is an error in eq (24), (1-exp(-2*a)) must be replaced by exp(-a) 
a = csi*h;
% % A = exp(-a).*(a.*PSI+(1+a).*PHI);
% % B = exp(-a).*((1-a).*PSI-a.*PHI);
A = (1-exp(-2.*a)).*(a.*PSI+(1+a).*PHI);
B = (1-exp(-2.*a)).*((1-a).*PSI-a.*PHI);
% *************************************************************************

% compute Uz and Ur, Fialko et al. (2001), eq. (12) and (13) **************
% NOTE there are two errors in eq (12) and (13)
% (1) 2*Uz and 2*Ur must be replaced by Uz and Ur
% (2) dcsi/sinh(csi*h) must be replaced by dcsi
Uz = zeros(size(r));                                                        % pre-allocate variable
Ur = zeros(size(r));                                                        % pre-allocate variable
for i=1:length(r)
    J0 = besselj(0,r(i)*csi);
    Uzi = J0.*(((1-2*nu)*B - csi*(z+h).*A).*sinh(csi*(z+h)) + ... 
                        (2*(1-nu)*A - csi*(z+h).*B).*cosh(csi*(z+h)));
    % % Uz(i) = wcsi*Uzi;  
    Uz(i) = 1/2.* wcsi*(Uzi.*1./sinh(a));  
    J1 = besselj(1,r(i)*csi);
    Uri = J1.*(((1-2*nu)*A + csi*(z+h).*B).*sinh(csi*(z+h)) + ... 
                        (2*(1-nu)*B + csi*(z+h).*A).*cosh(csi*(z+h)));
    % % Ur(i) = wcsi*Uri;    
    Ur(i) = 1/2.*wcsi*(Uri.*1./sinh(a));    
end                   
% *************************************************************************

% Deformation components **************************************************
u = rd*P_G*Ur.*x./r;
v = rd*P_G*Ur.*y./r;
w = -rd*P_G*Uz;
% *************************************************************************

% Volume change ***********************************************************
dV = -4*pi*(1-nu)*P_G*rd^3*(t*(wt'.*phi));
% *************************************************************************

function [pN, dpN] = legpol(x,N)
% legendre polynomial, Bonnet’s recursion formula
%
% Reference 
% http://en.wikipedia.org/wiki/Legendre_polynomials
%
% Note
% index n in wiki goes from 0 to N, in MATLAB j goes from 1 to N+1. To
% obtain the right coefficients we introduced j = n+1 -> n = j-1

P(1,:) = ones(size(x)); dP(1,:) = zeros(size(x));
P(2,:) = x; 
for j=2:N                                                               % loop up the recursion relation 
    P(j+1,:) = ((2*j-1)*x.*P(j,:) - (j-1)*P(j-1,:))/j;
     dP(j,:) = (j-1)*(x.*P(j,:) - P(j-1,:))./(x.^2-1);
end
pN =  P(N,:);
dpN = dP(N,:); 

function P = P(h,x)
% INPUT
%
% h     dimensionless source depth
% x     dummy variable
%
% OUTPUT
% P1-P4 are expressions from Appendix A of Fialko et al (2001). Used in the 
% definition of the functions T1 - T4 (see T.m)
%

P(1,:) = (12*h^2-x.^2)./(4*h^2+x.^2).^3;                                    % P1(x), pg 189
P(2,:) = log(4*h^2+x.^2) + (8*h^4+2*x.^2*h^2-x.^4)./(4*h^2+x.^2).^2;        % P2(x), pg 189
P(3,:) = 2*(8*h^4-2*x.^2*h^2+x.^4)./(4*h^2+x.^2).^3;                        % P3(x), pg 189
P(4,:) = (4*h^2-x.^2)./(4*h^2+x.^2).^2;                                     % P4(x), pg 189

function [x, w] = gauleg(x1,x2,N)
% function [x w] = gauleg(x1,x2,N)
% Given the upper and lower limits of integration x1 and x2, and given N,
% this routine returns arrays x(1:n) and w(1:n) of length N, containing the
% abscissas x and weights w of the Gaussian- legendre n-point quadrature
% formula. 
% *************************************************************************
% References.
% Loosely based on SUBROUTINE gauleg (Numerical Recipes in Fortran 77, 4.5) 
% *************************************************************************

z = zeros(1,N);         % pre-allocate variable
xm = 0.5*(x2+x1);       % mid-point
xl = 0.5*(x2-x1);       % half-interval

for n=1:N               % loop over the desidered roots
    z(n) = cos(pi*(n-0.25)/(N+0.5));          % approximation of the ith root of the Legendre's polynomials
    z1 = 100*z(n);
    while abs(z1-z(n)) > eps                   % Newton's method
        [pN, dpN] = legpol(z(n),N+1);          % compute the Legendre's polynomial and its derivative
        z1 = z(n);
        z(n) = z1 - pN/dpN;
    end
end
[~, dpN] = legpol(z,N+1);               % compute the derivative of the Legendre's polynomial
x(1:N) = xm - xl*z;                     % scale the root to the desidered interval
w(1:N) = 2*xl./((1-z.^2).*dpN.^2);      % compute the weights


function [phi, psi, t, w] = psi_phi(h)
% function [phi psi] = psi_phi(h,t)
% compute the function phi(t) and psi(t) using the Nystrom routine with the
% N-point Gauss-Legendre rule (Numerical Recipes in Fortran 77 (1992), 18.1
% Fredholm Equations of the Second Kind, p. 782.)
%
% INPUT
% h     dimensionless source depth
% t     dummy variable vector (value between 0 and 1; length(t) = M)
%
% OUTPUT
% phi and psi are expressions from Appendix A of Fialko et al (2001), 
% equation (A1). phi and psi are used in the definition of the functions 
% PHI and PSI. See also A_B.m and Fialko et al (2001), eq. (26), p. 184.
% t and w are the abscissas and weights from the Gauss-Legendre quadrature
% rule
%

[t, w] = gauleg(0,1,41);                                                     % abscissas and weights from the Gauss-Legendre quadrature rule

% Solution at the quadrature points r(j) 
g = -2*t/pi;                                                                % see eq. (A1) in Fialko et al. (2001)
d = [g zeros(size(g))]';                                                     

[T1, T2, T3, T4] = T(h,t,t);                                                   % Fredholm's integration kernels (Falko et al., 2001, eq. 27)
T1tilde = zeros(size(T1)); T2tilde = zeros(size(T1));                       % pre-allocate variable
T3tilde = zeros(size(T1)); T4tilde = zeros(size(T1));

N = length(t); 
for j= 1:N                                                                  % multiply the integration kernels by the quadrature weights
    T1tilde(:,j) = w(j)*T1(:,j);                                            
    T2tilde(:,j) = w(j)*T2(:,j);
    T3tilde(:,j) = w(j)*T3(:,j);
    T4tilde(:,j) = w(j)*T4(:,j);
end

Ktilde = [T1tilde, T3tilde; T4tilde, T2tilde];                              % see eq (18.1.5) in Press et al (1992)
y = (eye(2*N,2*N)-(2/pi)*Ktilde)\d;                                         % solution of linear system, (18.1.6) in Press et al (1992)
phi = y(1:N);                                                               % phi at the quadrature points t(j)
psi = y(N+1:2*N);                                                           % psi at the quadrature points t(j)


function [T1, T2, T3, T4] = T(h,t,r)
% INPUT
%
% h     dimensionless source depth
% r     dummy variable (used to integrate along the sill dimenionless radius) 
% t     dummy variable
%
% OUTPUT
% T1-T4 are expressions from Appendix A of Fialko et al (2001), equations 
% (A2) - (A5). T1-T4 are used in the definition of the functions phi and 
% psi (see phi_psi.m)

M = length(t); N = length(r); 
T1 = zeros(M,N); T2 = T1; T3 = T1;                                          % pre-allocate variables

for i=1:M
    Pm = P(h,t(i)-r);                                                       % define functions P1-P4
    Pp = P(h,t(i)+r);                                                       % define functions P1-P4
    T1(i,:) = 4*h^3*(Pm(1,:)-Pp(1,:));                                      % equation (A2)
    T2(i,:) = (h./(t(i).*r)).*(Pm(2,:)-Pp(2,:)) +h*(Pm(3,:)+Pp(3,:));       % equation (A3)
    T3(i,:) = (h^2./r).*(Pm(4,:)-Pp(4,:)...
                              -2*r.*((t(i)-r).*Pm(1,:)+(t(i)+r).*Pp(1,:))); % equation (A4)
end
T4 = T3';                                                               % equation (A5)
