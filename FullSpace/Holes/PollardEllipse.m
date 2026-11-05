function [S, U, C] = PollardEllipse(x, y, p, s1, s3, beta, a, b, G, nu, varargin)
%POLLARD_ELLIPSE  Stress and displacement around a pressurized elliptical hole.
%
%   [S,U,C] = POLLARD_ELLIPSE(X,Y,P,S1,S3,BETA,A,B,G,NU,...)
%
%   Plane-elasticity solution for an elliptical hole in an infinite,
%   homogeneous, isotropic elastic sheet, with a uniform normal stress P on
%   the hole wall and principal stresses S1, S3 at infinity.  Follows
%   Pollard, D.D., 1973, "Equations for stress and displacement fields
%   around pressurized elliptical holes in elastic solids", Mathematical
%   Geology 5(1), 11-25 (complex potentials of Kolosoff, 1914; eqs 34-42).
%
%   The stresses/displacements are evaluated here directly from the complex
%   potentials rather than from the expanded real forms (eqs 43-57); the
%   result is identical but shorter and better behaved numerically.
%
%   INPUTS
%     X,Y    evaluation points; equal-size arrays (or scalars, expanded).
%            Origin at the hole centre, +x along the long axis.
%            Points inside the hole return NaN.
%     P      uniform normal stress on the hole wall (sigma_xixi = P there).
%     S1,S3  principal stresses at infinity.  S3 acts at angle BETA to the
%            +x axis, S1 acts at BETA + pi/2  (eqs 29-31).
%     BETA   orientation of S3 from +x, radians by default.
%     A,B    semi-major (x) and semi-minor (y) axes of the hole, i.e. l/2
%            and t/2 in Pollard's notation.  B = 0 gives a flat crack.
%            Requires A > B >= 0.
%     G      shear modulus.
%     NU     Poisson's ratio.
%
%   OPTIONS  (name/value pairs)
%     'plane'       'strain' (default) | 'stress'
%     'convention'  'tension'     (default) tension-positive, the convention
%                                 the published equations are written in: a
%                                 fluid pressure that pushes the walls apart
%                                 is a NEGATIVE P here.
%                   'compression' compression-positive (usual geological
%                                 convention): input and output stresses
%                                 both change sign, so P > 0 is an outward
%                                 fluid pressure that opens the hole and
%                                 S1 >= S3 are compressions.  Displacements
%                                 are physical vectors either way.
%     'angleunits'  'rad' (default) | 'deg'   (applies to BETA and to the
%                                 output angles ALPHA and GAMMA)
%     'remotestrain' 'include' (default) | 'remove'
%                   'remove' subtracts the asymptotic linear displacement
%                   field (the uniform remote strain plus the accompanying
%                   rigid rotation), leaving only the hole-induced
%                   perturbation, which decays to zero at infinity.  Use
%                   this for opening/closing of a crack or dike under
%                   remote load, where the uniform far-field stretch of the
%                   whole sheet is not wanted.  Stresses are unaffected;
%                   with no remote load (S1 = S3 = 0) the displacements are
%                   unaffected too, since the asymptotic field is then zero.
%
%   OUTPUTS
%     S  struct of stresses, each the size of X:
%          .xx .yy .xy        Cartesian components
%          .xixi .etaeta .xieta   elliptical components (eqs 21-23)
%          .p1 .p3            in-plane principal stresses (p1 >= p3 for
%                             'tension'; p1 is the most compressive for
%                             'compression')
%          .tmax              maximum in-plane shear stress
%          .alpha             angle from +x to the p1 direction (CCW)
%     U  struct of displacements:
%          .x .y              Cartesian components (eqs 37-38)
%          .xi .eta           elliptical components (eqs 24-25)
%          .mag               magnitude
%          .gamma             angle from +x to the displacement vector
%     C  struct of auxiliary quantities:
%          .xi .eta           elliptical coordinates of each point (eqs 7-8)
%          .f                 focal half-distance, sqrt(A^2-B^2)
%          .xi0               ellipse defining the hole, atanh(B/A) (eq 27)
%          .kappa             3-4*NU (plane strain) or (3-NU)/(1+NU)
%
%   EXAMPLES
%     % Hoop stress at the tip of a 4:1 hole under remote tension Syy = 1
%     % (Inglis: 1 + 2a/b = 9)
%     [S,~,~] = pollard_ellipse(2, 0, 0, 1, 0, 0, 2, 0.5, 1, 0.25);
%     S.etaeta                                  % -> 9
%
%     % Opening of a fluid-pressurized crack, geological sign convention
%     a = 1; b = 1e-3; th = linspace(0,pi,200);
%     [~,U] = pollard_ellipse(a*cos(th), b*sin(th), 10e6, 0, 0, 0, ...
%                             a, b, 30e9, 0.25, 'convention','compression');
%     plot(a*cos(th), U.y)                      % half-opening profile
%
%     % Field on a grid
%     [X,Y] = meshgrid(linspace(-4,4,201), linspace(-3,3,151));
%     S = pollard_ellipse(X, Y, 5, 2, -1, pi/6, 2, 0.6, 30e3, 0.25);
%     pcolor(X,Y,S.yy); shading interp; axis equal; colorbar

% ---------------------------------------------------------------- options
opt.plane      = 'strain';
opt.convention = 'tension';
opt.angleunits = 'rad';
opt.remotestrain = 'include';
if mod(numel(varargin),2) ~= 0
    error('pollard_ellipse:options','Options must be name/value pairs.');
end
for k = 1:2:numel(varargin)
    name = lower(varargin{k});
    if ~isfield(opt,name)
        error('pollard_ellipse:options','Unknown option ''%s''.',varargin{k});
    end
    opt.(name) = lower(varargin{k+1});
end
switch opt.plane
    case 'strain', kappa = 3 - 4*nu;
    case 'stress', kappa = (3 - nu)/(1 + nu);
    otherwise, error('pollard_ellipse:plane','''plane'' must be ''strain'' or ''stress''.');
end
switch opt.convention
    case 'tension',     sgn =  1;
    case 'compression', sgn = -1;
    otherwise, error('pollard_ellipse:convention','''convention'' must be ''tension'' or ''compression''.');
end
switch opt.remotestrain
    case 'include', drop_remote = false;
    case 'remove',  drop_remote = true;
    otherwise, error('pollard_ellipse:remotestrain','''remotestrain'' must be ''include'' or ''remove''.');
end
switch opt.angleunits
    case 'rad', aunit = 1;
    case 'deg', aunit = pi/180;
    otherwise, error('pollard_ellipse:angleunits','''angleunits'' must be ''rad'' or ''deg''.');
end

% ------------------------------------------------------------ input check
if ~(isscalar(a) && isscalar(b) && a > 0 && b >= 0 && b < a)
    error('pollard_ellipse:axes', ...
        ['Require A > B >= 0 (semi-major along x, semi-minor along y). ' ...
         'For a circular hole use B = A*(1-1e-6) or the Kirsch solution.']);
end
if ~(isscalar(G) && G > 0), error('pollard_ellipse:G','G must be a positive scalar.'); end
if ~(isscalar(nu) && nu > -1 && nu < 0.5), error('pollard_ellipse:nu','NU must be in (-1,0.5).'); end
if isscalar(x) && ~isscalar(y), x = x*ones(size(y)); end
if isscalar(y) && ~isscalar(x), y = y*ones(size(x)); end
if ~isequal(size(x),size(y)), error('pollard_ellipse:size','X and Y must be the same size.'); end

% stresses to the tension-positive convention the equations are written in
p    = sgn*p;
s1   = sgn*s1;
s3   = sgn*s3;
beta = beta*aunit;

% ---------------------------------------------- elliptical coordinates (7,8)
f    = sqrt(a^2 - b^2);            % focal half-distance
xi0  = atanh(b/a);                 % ellipse bounding the hole, eq (27)
z    = x + 1i*y;
zeta = acosh(z/f);                 % z = f*cosh(zeta), eq (17)
zeta(real(zeta) < 0) = -zeta(real(zeta) < 0);   % keep xi >= 0 (cosh is even)
xi   = real(zeta);
eta  = imag(zeta);
sh   = sinh(zeta);
ch   = cosh(zeta);

% ------------------------------------------- complex potentials, eqs (39-42)
E     = exp(2*xi0);
c2xi0 = cosh(2*xi0);
psi = zeros(size(z)); psi_z = psi; psi_zz = psi; chi_z = psi; chi_zz = psi;
Gam = 0; Gamp = 0;   % coefficients of the asymptotic forms psi -> Gam*z,
                     % dchi/dz -> Gamp*z  (used only by 'remotestrain')

% solution (1): uniaxial s3-p at angle beta;  solution (2): s1-p at beta+pi/2
Amp = [s3 - p, s1 - p];
Ang = [beta,   beta + pi/2];
for k = 1:2
    A = Amp(k);  B = Ang(k);
    if A == 0, continue; end
    Gam  = Gam  + (A/4)*(1 - 1i*E*sin(2*B));
    Gamp = Gamp - (A/2)*exp(-2i*B);
    c2b  = cos(2*B);
    e2ib = exp(2i*B);
    k1   = A*f/4;
    P1   = k1*(E*c2b*ch + (1 - E*e2ib).*sh);          % psi_k          (39)
    P1z  = k1*(E*c2b*sh + (1 - E*e2ib).*ch);          % d psi_k / dzeta
    k2   = -A*f^2/4;
    arg  = 2*(zeta - xi0 - 1i*B);
    X1z  = k2*((c2xi0 - c2b) + E*sinh(arg));          % d chi_k / dzeta (40)
    X1zz = k2*(2*E*cosh(arg));                        % d2chi_k/ dzeta2
    psi    = psi    + P1;
    psi_z  = psi_z  + P1z;
    psi_zz = psi_zz + P1;                             % d2psi/dzeta2 = psi
    chi_z  = chi_z  + X1z;
    chi_zz = chi_zz + X1zz;
end
% solution (3): hydrostatic p over the whole sheet, eqs (41-42)
psi    = psi    + p*f*ch/2;
psi_z  = psi_z  + p*f*sh/2;
psi_zz = psi_zz + p*f*ch/2;
Gam    = Gam    + p/2;

% ------------------------------------- derivatives w.r.t. z  (dz/dzeta = f*sh)
d1psi = psi_z./(f*sh);
d2psi = (psi_zz.*sh - psi_z.*ch)./(f^2*sh.^3);
d1chi = chi_z./(f*sh);
d2chi = (chi_zz.*sh - chi_z.*ch)./(f^2*sh.^3);

% ------------------------------------------- Cartesian stresses, eqs (34-36)
T   = conj(z).*d2psi + d2chi;
sxx = real(2*d1psi - T);
syy = real(2*d1psi + T);
sxy = imag(T);

% ---------------------------------------- Cartesian displacements, eqs (37-38)
uu = (kappa*psi - z.*conj(d1psi) - conj(d1chi))/(2*G);
if drop_remote
    % asymptotic (uniform strain + rigid rotation) part of the field
    uu = uu - ((kappa*Gam - conj(Gam))*z - conj(Gamp)*conj(z))/(2*G);
end
ux = real(uu);
uy = imag(uu);

% ------------------------ rotation into elliptical components, eqs (9-11,19-25)
e2id = sh./conj(sh);                       % exp(2*i*delta), eq (19)
W    = e2id.*((syy - sxx) + 2i*sxy);       % eq (10)
sxixi   = ((sxx + syy) - real(W))/2;
setaeta = ((sxx + syy) + real(W))/2;
sxieta  = imag(W)/2;
eid  = sh./abs(sh);                        % exp(i*delta)
uc   = conj(eid).*uu;                      % eq (11)
uxi  = real(uc);
ueta = imag(uc);

% ------------------------------------------------ blank out points in the hole
inside = xi < xi0 - 1e-9*max(1,xi0);
sxx(inside)     = NaN;  syy(inside)     = NaN;  sxy(inside)    = NaN;
sxixi(inside)   = NaN;  setaeta(inside) = NaN;  sxieta(inside) = NaN;
ux(inside)      = NaN;  uy(inside)      = NaN;
uxi(inside)     = NaN;  ueta(inside)    = NaN;

% ------------------------------------------------------- back to user's signs
sxx = sgn*sxx; syy = sgn*syy; sxy = sgn*sxy;
sxixi = sgn*sxixi; setaeta = sgn*setaeta; sxieta = sgn*sxieta;

% --------------------------------------- principal stresses and angles [38-41]
mean_s = (sxx + syy)/2;
tmax   = sqrt(((sxx - syy)/2).^2 + sxy.^2);
p1     = mean_s + tmax;
p3     = mean_s - tmax;
alpha  = 0.5*atan2(2*sxy, sxx - syy);
alpha(tmax < eps*max(1,max(abs(mean_s(:))))) = NaN;   % hydrostatic: undefined
umag   = hypot(ux, uy);
gamma  = atan2(uy, ux);
gamma(umag == 0) = NaN;

% ----------------------------------------------------------------- pack output
S = struct('xx',sxx,'yy',syy,'xy',sxy, ...
           'xixi',sxixi,'etaeta',setaeta,'xieta',sxieta, ...
           'p1',p1,'p3',p3,'tmax',tmax,'alpha',alpha/aunit);
U = struct('x',ux,'y',uy,'xi',uxi,'eta',ueta,'mag',umag,'gamma',gamma/aunit);
C = struct('xi',xi,'eta',eta,'f',f,'xi0',xi0,'kappa',kappa);
end
