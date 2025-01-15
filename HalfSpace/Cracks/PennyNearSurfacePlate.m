function [uz,ur] = PennyNearSurfacePlate(d,a,nu,E,p,r)
    %Asymptotic solution for a penny-shaped
    %near-surface hydraulic fracture
    %a/d > 5 -works well for this...
%		  r: radial distance of observation,
%		  d: depth of the center of the source from the surface,
%		  a: radius of the source with the hydrostatic pressure,
%		  p: change of the hydrostatic pressure in the crack.
%		  E: Young's modulus,
%		 nu: Poisson's ratio (default is 0.25 for isotropic medium).    



    H=d;
    R=a;
    rho=r./R;
    D=(E/(1-nu^2)*H^3)/12;%Eq.1
    delta=0.620;%eq.15
    eta=(4.*delta)./(3+2.*nu).*H./R;%Eq.21
    uz=(R.^4.*p)./(64.*D).*((1-rho.^2).^2+4.*eta.*(1-rho.^2));%Eq.22

    % syms R p D rho eta r
    % rho=r./R;
    % uz2=(R.^4.*p)./(64.*D).*((1-rho.^2).^2+4.*eta.*(1-rho.^2));
    % diff(uz2)/diff(r)
    dwdr=(R^4.*p.*((4.*r.*(r.^2./R.^2 - 1))./R.^2 - (8.*eta.*r)./R.^2))./(64.*D);
    ur=-(dwdr/2*d);



end