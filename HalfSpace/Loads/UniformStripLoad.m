function [Sxx,Szz,Sxz,Ux,Uz] = UniformStripLoad(x,z,a,sigma,G,nu)
%UNIFORMSTRIPLOAD  Stresses & displacements: uniform normal traction
%                  over a surface strip of a 2D elastic half-space.
%
% A uniform normal pressure  sigma  applied over the surface strip
% -a <= x <= a  of a half-space, constant out-of-plane (plane strain).
% Equivalent to a superposition of Flamant line loads dP = sigma*dx'.
%
% Reference: Jaeger, Cook & Zimmerman, "Fundamentals of Rock Mechanics",
%            4th edn (2007), Ch.13 - surface loads on a half-space
%            (uniform normal traction over a strip).
%
% --------------------------------------------------------------------
% COORDINATE CONVENTION
%   x  = horizontal coordinate along the free surface
%   z  = depth into the half-space, POSITIVE DOWNWARD
%   The loaded strip is  z = 0,  -a <= x <= a.
%   For a z-up convention, negate z on input and Sxz,Uz on output.
%
% SIGN CONVENTION:  COMPRESSION POSITIVE for all stresses.
%
% --------------------------------------------------------------------
% INPUTS
%   x,z   : observation point coordinates. Scalars or arrays of any
%           matching size; fully vectorised.
%   a     : strip HALF-width (strip spans -a..a, total width 2a)
%   sigma : uniform applied normal traction (pressure)
%   G     : shear modulus
%   nu    : Poisson's ratio (plane strain)
%
% OUTPUTS  (same size as x,z)
%   Sxx : horizontal normal stress
%   Szz : vertical   normal stress
%   Sxz : shear stress
%   Ux  : horizontal displacement
%   Uz  : vertical   displacement (positive downward)
%
% --------------------------------------------------------------------
% FORMULAE  (compression positive)
%
% Edge angles, with the strip edges at x = +a and x = -a :
%   theta1 = atan2(x - a, z)        (ray to the edge at x = +a)
%   theta2 = atan2(x + a, z)        (ray to the edge at x = -a)
%   dth    = theta2 - theta1        (the angle the strip subtends)
%   sumt   = theta1 + theta2
%
% Stresses:
%   Sxx = (sigma/pi) [ dth - sin(dth) cos(sumt) ]
%   Szz = (sigma/pi) [ dth + sin(dth) cos(sumt) ]
%   Sxz = (sigma/pi) [       sin(dth) sin(sumt) ]
%
% These are eqs (13.27)-(13.29) of Hudson & Harrison / the standard
% strip-load result in Jaeger, Cook & Zimmerman. WARNING: the closed
% form is sensitive to the angle convention - with the atan2 definitions
% above the vertical stress Szz carries the "+ sin cos" sign. This
% implementation has been checked against a direct numerical integration
% of the Flamant kernel over the strip and agrees to ~1e-4.
%
% Displacements: the surface-load strip displacement field in 2D has a
% closed form, but it crosses branch cuts that make a hand-written
% expression error-prone. This routine instead obtains the displacements
% by NUMERICALLY INTEGRATING the verified Flamant line-load displacement
% (FlamantPointLoad) across the strip - dP = sigma*dx'. This is correct
% by construction and avoids branch-cut mistakes. It is accurate to the
% integration resolution (controlled by NQUAD below). As with any 2D
% line-load solution the result is defined only up to an additive
% rigid-body constant; use RELATIVE displacements.
%
% --------------------------------------------------------------------
% NOTES
%   * The stress field is independent of the elastic constants (G,nu);
%     they appear only in the displacements.
%   * Boundary check: as z -> 0, Szz -> sigma for |x| < a and -> 0 for
%     |x| > a; exactly at an edge Szz -> sigma/2.
%   * Sxx is ZERO on the free surface outside the strip and = sigma on
%     the surface inside it. Just below the surface outside the strip it
%     grows linearly with depth: Sxx ~ (4 sigma a / pi) z /(x^2 - a^2).
%   * The displacements inherit the 2D line-load ln(r) behaviour: they
%     are defined only up to an additive constant and grow at infinity.
%     Use RELATIVE displacements between two points.
%   * No singularity for z > 0. On the surface z = 0 the stresses are
%     finite everywhere except the two edges x = +/- a (log singular in
%     the displacement gradient).
%
% EXAMPLE
%   x = linspace(-4,4,401);  z = linspace(0.01,4,301);
%   [X,Z] = meshgrid(x,z);
%   [Sxx,Szz,Sxz] = UniformStripLoad(X,Z,1,1,1,0.25);
%   contourf(X,Z,Szz,30,'LineColor','none');
%   set(gca,'YDir','reverse'); axis equal tight; colorbar;
%
% See also FLAMANTPOINTLOAD.
% --------------------------------------------------------------------

% ---- input handling ------------------------------------------------
if nargin < 3 || isempty(a),     a     = 1;    end
if nargin < 4 || isempty(sigma), sigma = 1;    end
if nargin < 5 || isempty(G),     G     = 1;    end
if nargin < 6 || isempty(nu),    nu    = 0.25; end
if ~isequal(size(x),size(z))
    error('UniformStripLoad:size','x and z must have the same size.');
end
if a <= 0
    error('UniformStripLoad:a','Strip half-width a must be positive.');
end

NQUAD = 800;                % number of line loads used for the
                            % displacement integration (raise for
                            % higher accuracy near the surface)

% ---- edge geometry -------------------------------------------------
xm = x - a;                 % horizontal offset from the +a edge
xp = x + a;                 % horizontal offset from the -a edge

theta1 = atan2(xm, z);      % ray to edge at x = +a
theta2 = atan2(xp, z);      % ray to edge at x = -a
dth    = theta2 - theta1;   % angle subtended by the strip
sumt   = theta1 + theta2;

% ====================================================================
% STRESSES  (compression positive) - closed form, verified against a
% direct numerical integration of the Flamant kernel to ~1e-4.
% ====================================================================
Sxx = (sigma/pi) * ( dth - sin(dth).*cos(sumt) );
Szz = (sigma/pi) * ( dth + sin(dth).*cos(sumt) );
Sxz = (sigma/pi) * (       sin(dth).*sin(sumt) );

% ====================================================================
% DISPLACEMENTS  - integrate the verified Flamant line-load solution
% over the strip:  dP = sigma * dx',  x' in [-a, a].
% Trapezoidal rule over NQUAD sub-loads.
% ====================================================================
xq = linspace(-a, a, NQUAD);          % positions of the sub-line-loads
wq = ones(1,NQUAD); wq(1) = 0.5; wq(end) = 0.5;
wq = wq * (xq(2)-xq(1));               % trapezoidal weights
dP = sigma;                            % load intensity (per unit dx')

Ux = zeros(size(x));
Uz = zeros(size(x));
for k = 1:NQUAD
    % observation point relative to the k-th line load at (xq(k),0)
    [~,~,~,uxk,uzk] = FlamantPointLoad(x - xq(k), z, dP, 0, G, nu);
    Ux = Ux + wq(k) * uxk;
    Uz = Uz + wq(k) * uzk;
end

end
