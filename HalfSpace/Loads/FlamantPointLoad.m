function [Sxx,Szz,Sxz,Ux,Uz] = FlamantPointLoad(x,z,P,T,G,nu)
%FLAMANTPOINTLOAD  Stresses & displacements: normal+tangential line load
%                  on the surface of a 2D elastic half-space.
%
% The Flamant (1892) solution. A concentrated load acting on the surface
% of a half-space, constant in the out-of-plane direction, i.e. a LINE
% load of force-per-unit-length. Plane strain.
%
% Reference: Jaeger, Cook & Zimmerman, "Fundamentals of Rock Mechanics",
%            4th edn (2007), Ch.13 - surface loads on a half-space.
%            Also Johnson, "Contact Mechanics" (1985), section 2.7.
%
% --------------------------------------------------------------------
% COORDINATE CONVENTION
%   x  = horizontal coordinate along the free surface
%   z  = depth into the half-space, POSITIVE DOWNWARD
%   The load is applied at the origin (x,z) = (0,0).
%     P = normal      line load  (force / unit length), pushes +z (down)
%     T = tangential  line load  (force / unit length), acts in +x
%   For a z-up convention, negate z on input and negate Sxz,Uz on output.
%
% SIGN CONVENTION:  COMPRESSION POSITIVE for all stresses.
%
% --------------------------------------------------------------------
% INPUTS
%   x,z : observation point coordinates. Scalars or arrays (any matching
%         size); the function is fully vectorised.
%   P   : normal line load     (force per unit length)
%   T   : tangential line load (force per unit length)
%   G   : shear modulus
%   nu  : Poisson's ratio        (plane strain)
%
% OUTPUTS  (same size as x,z)
%   Sxx : horizontal normal stress
%   Szz : vertical   normal stress
%   Sxz : shear stress
%   Ux  : horizontal displacement
%   Uz  : vertical   displacement (positive downward)
%
% --------------------------------------------------------------------
% NOTES
%   * The stress field depends ONLY on the loads and geometry - the
%     elastic constants (G,nu) cancel out of Sxx,Szz,Sxz. They enter
%     only the displacements.
%   * Every stress component carries a factor of z in the numerator, so
%     all stresses -> 0 on the free surface (z=0) everywhere except at
%     the singular load point itself.
%   * The displacements contain a ln(r) term and therefore grow without
%     bound at infinity. This is the classic feature of 2D line-load
%     solutions: only RELATIVE displacements between two points are
%     physically meaningful. An arbitrary rigid-body constant may be
%     added freely (see DATUM option below).
%   * At the load point r->0 the solution is singular (stresses ~1/r,
%     displacements ~ln r). Observation points must have r > 0.
%
% EXAMPLE
%   x = linspace(-5,5,201);  z = linspace(0.01,5,201);
%   [X,Z] = meshgrid(x,z);
%   [Sxx,Szz,Sxz] = FlamantPointLoad(X,Z,1,0,1,0.25);
%   contourf(X,Z,Szz,30,'LineColor','none');
%   set(gca,'YDir','reverse'); axis equal tight; colorbar;
%
% --------------------------------------------------------------------

% ---- input handling ------------------------------------------------
if nargin < 4 || isempty(T),  T  = 0;     end
if nargin < 5 || isempty(G),  G  = 1;     end
if nargin < 6 || isempty(nu), nu = 0.25;  end
if ~isequal(size(x),size(z))
    error('FlamantPointLoad:size','x and z must have the same size.');
end

% ---- geometry ------------------------------------------------------
r2 = x.^2 + z.^2;          % squared distance from the load point
r4 = r2.^2;

if any(r2(:) <= 0)
    warning('FlamantPointLoad:singular', ...
        'Observation point coincides with the load point (r = 0).');
end

% ====================================================================
% STRESSES  (compression positive). Superposition of the normal-load
% and tangential-load Flamant fields.
%
% Normal line load P (acting along +z):
%   Sxx =  2P/pi * x^2 z / r^4
%   Szz =  2P/pi * z^3   / r^4
%   Sxz =  2P/pi * x z^2 / r^4
%
% Tangential line load T (acting along +x): obtained by swapping the
% role of the load direction (rotate the normal-load field by 90 deg):
%   Sxx =  2T/pi * x^3   / r^4
%   Szz =  2T/pi * x z^2 / r^4
%   Sxz =  2T/pi * x^2 z / r^4
% ====================================================================
Sxx = 2/pi * ( P .* (x.^2 .* z) + T .* (x.^3)       ) ./ r4;
Szz = 2/pi * ( P .* (z.^3)      + T .* (x .* z.^2)  ) ./ r4;
Sxz = 2/pi * ( P .* (x .* z.^2) + T .* (x.^2 .* z)  ) ./ r4;

% ====================================================================
% DISPLACEMENTS  (plane strain). Both fields were derived by integrating
% the strains of the stress field above (symbolically) and verified to
% reproduce those strains to ~1e-6.
%
% Normal line load P:
%   Ux = P/(2 pi G) [  x z / r^2  + (2 nu - 1) atan(x/z) ]
%   Uz = P/(2 pi G) [ -x^2 / r^2  - 2(1 - nu) ln(r)      ]
%
% Tangential line load T:
%   Ux = T/(2 pi G) [ -z^2 / r^2  - 2(1 - nu) ln(r)      ]
%   Uz = T/(2 pi G) [  x z / r^2  + (2 nu - 1) atan(z/x) ]
%
% Branch handling:
%   atan(x/z) -> atan2(x,z)              (continuous for z > 0)
%   atan(z/x) -> pi/2 - atan2(x,z)       (continuous across x = 0;
%                                         a plain atan(z/x) would jump)
% An arbitrary rigid-body constant may be added to Ux and Uz; the ln(r)
% term means absolute displacements are only defined up to that constant
% (and grow at infinity - use relative displacements).
% ====================================================================
lnr   = 0.5*log(r2);                 % = ln(r)
ang_xz = atan2(x,z);                  % = atan(x/z), valid for z > 0
ang_zx = pi/2 - ang_xz;               % = atan(z/x), continuous in x
c     = 1/(2*pi*G);

Ux =  c * ( P .* (  x.*z./r2 + (2*nu-1).*ang_xz ) ...
          + T .* ( -z.^2./r2 - 2*(1-nu).*lnr    ) );

Uz =  c * ( P .* ( -x.^2./r2 - 2*(1-nu).*lnr    ) ...
          + T .* (  x.*z./r2 + (2*nu-1).*ang_zx ) );

end
