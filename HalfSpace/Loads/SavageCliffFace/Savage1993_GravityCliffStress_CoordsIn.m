function [Sxx,Syy,Sxy,x,y] = Savage1993_GravityCliffStress_CoordsIn(rg,nu,mu,b,x,y)
% Savage1993_GravityCliffStress_CoordsIn: Returns Cartesian stresses beneath
%              a vertical cliff (bench of height b) loaded by gravity, at
%              real x,y locations. Same solution as
%              Savage1993_GravityCliffStress; here the conformal map of
%              Savage (1993) is inverted numerically (Newton) to find u,v.
%
%              Geometry: top of the cliff y = 0 for x >= 0, cliff face x = 0
%              for -b <= y <= 0, lower ground y = -b for x <= 0. Points in
%              the air return NaN.
%
%              Savage, W.Z., 1993. Gravity-induced stresses near a vertical
%              cliff. Int. J. Rock Mech. Min. Sci. & Geomech. Abstr. 30(4),
%              325-330.
%
% usage #1:
% [Sxx,Syy,Sxy,x,y] = Savage1993_GravityCliffStress_CoordsIn(rg,nu,mu,b,x,y)
%
% Arguments: (input)
%       rg    - Weight, the density of material times the acceleration due to
%               gravity.
%
%       nu    - Poisson's Ratio.
%
%       mu    - Shear Modulus (not used; stresses do not depend on it).
%
%       b     - Bench (cliff) height.
%
%       x,y   - Observation points, actual x and y locations.
%
% Arguments: (output)
%       Sxx,Syy,Sxy - 2D stress tensor components, total stress (gravity),
%               tension positive. NaN in the air (and exactly at the base
%               corner, where the stress is singular).
%
%       x,y   - same as inputs
%
% Example usage:
%
%  nu = 0.25;
%  mu = 500;
%  rg = 2700*9.81;
%  b  = 100;
%  [x,y] = meshgrid(linspace(-100,150,251),linspace(-300,0,301));
%
%  [Sxx,Syy,Sxy] =...
%  Savage1993_GravityCliffStress_CoordsIn(rg,nu,mu,b,x,y);
%
%  DrawScatterPlots2d( x,y,[],Sxx,Syy,Sxy )
%
%  Needs Savage1993_GravityCliffStress.m
%  Written with Claude (Anthropic), Nov 2026, for Tim Davis' ElasticSolutions.

sz  = size(x);
z   = complex(x(:), y(:))/b;                       % work with b = 1
tol = 1e-12;
inrock = ((imag(z) <= tol) & (real(z) >= -tol)) | (imag(z) <= -1 + tol);
inrock = inrock & abs(z + 1i) > 1e-12;            % base corner is singular -> NaN
Sxx = NaN(numel(z),1);  Syy = Sxx;  Sxy = Sxx;

if any(inrock)
    w = invertmap(z(inrock));
    ok = isfinite(w);
    idx = find(inrock);
    [a, c, e] = Savage1993_GravityCliffStress(rg, nu, mu, b, real(w(ok)), imag(w(ok)));
    Sxx(idx(ok)) = a;  Syy(idx(ok)) = c;  Sxy(idx(ok)) = e;
end
Sxx = reshape(Sxx, sz);  Syy = reshape(Syy, sz);  Sxy = reshape(Sxy, sz);
end

% =========================================================================
function w = invertmap(z)
% Solve f(w) = z (b = 1), v <= 0, by damped Newton. Starting guesses come
% from the local forms of the map, and the best of them is used for each point:
%   far field   z ~ (w + log(2w))/pi
%   top corner  z ~ (2 sqrt(2)/pi) sqrt(w-1)
%   base corner z + i ~ (sqrt(2)/(3 pi)) i (w+1)^(3/2)
n  = numel(z);
zp = z + 1i;
% far field: fixed-point iterations on w = pi z - log(2w)
wf = pi*z;
for k = 1:6
    wf = pi*z - log(2*wf);
    wf = complex(real(wf), min(imag(wf), -1e-12));
end
% top corner
wt = 1 + (pi*z/(2*sqrt(2))).^2;
% base corner: angle of z+i taken in (-5pi/4, 3pi/4], body occupies [-pi, pi/2]
th = angle(zp);  th(th > 0.75*pi) = th(th > 0.75*pi) - 2*pi;
wb = -1 + (3*pi*sqrt(2)/2*abs(zp)).^(2/3).*exp(1i*(2/3)*(th - pi/2));
W0 = [wf, wt, wb];
W0 = complex(real(W0), min(imag(W0), -1e-15));
R0 = abs(fmap(W0) - z);
[~, best] = min(R0, [], 2);
w  = W0(sub2ind(size(W0), (1:n).', best));

[fz, df] = fmap(w);
res = abs(fz - z);
tolr = 1e-13*max(1, abs(z)) + 1e-14*abs(df);
act = find(res > tolr);
for it = 1:80
    if isempty(act), break; end
    dw  = (z(act) - fz(act))./df(act);
    lam = ones(numel(act),1);
    wn = w(act); fn = fz(act); dn = df(act); rn = res(act);
    todo = (1:numel(act)).';
    for k = 1:30                                   % halve the step until the residual drops
        wt = w(act(todo)) + lam(todo).*dw(todo);
        wt = complex(real(wt), min(imag(wt), -1e-15));
        [ft, dt] = fmap(wt);
        rt = abs(ft - z(act(todo)));
        ok = rt < res(act(todo));
        wn(todo(ok)) = wt(ok); fn(todo(ok)) = ft(ok); dn(todo(ok)) = dt(ok); rn(todo(ok)) = rt(ok);
        todo = todo(~ok);
        if isempty(todo), break; end
        lam(todo) = lam(todo)/2;
    end
    stalled = false(numel(act),1);  stalled(todo) = true;
    w(act) = wn; fz(act) = fn; df(act) = dn; res(act) = rn;
    tolr(act) = 1e-13*max(1, abs(z(act))) + 1e-14*abs(dn);
    act = act(res(act) > tolr(act) & ~stalled);
end
w(res > max(1e-8*max(1, abs(z)), 10*tolr)) = NaN;   % genuine failures only
end

function [z, g] = fmap(w)
% Savage's map with b = 1 (accurate near both corners) and its derivative
w  = complex(real(w), min(imag(w), -1e-15));
r2 = sqrt(2);
sm = sqrt(w - 1);  sp = sqrt(w + 1);
zb = (sm.*sp + 2*(log1pc((w + 1)./((sm - 1i*r2).*(-1i*r2))) + log1pc(sp./sm)))/pi - 1i;
zt = (sm.*sp + 2*(log1pc((w - 1)./((sp + r2)*r2)) + log1pc(sm./sp)))/pi;
z  = zb;  k = abs(w - 1) < abs(w + 1);  z(k) = zt(k);
g  = (sp./sm)/pi;
end

function y = log1pc(x)
u = 1 + x;  y = log(u);  k = (u ~= 1);
y(k) = y(k).*x(k)./(u(k) - 1);  y(~k) = x(~k);
end
