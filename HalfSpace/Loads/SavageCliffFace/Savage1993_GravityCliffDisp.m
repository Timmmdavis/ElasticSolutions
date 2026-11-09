function [Ux,Uy,X,Y] = Savage1993_GravityCliffDisp(rg,nu,mu,b,u,v,xref)
% Savage1993_GravityCliffDisp: Returns the displacements caused by
%              excavating a vertical cliff (bench of height b) in a body
%              loaded by gravity. Points are given in the conformally mapped
%              (u,v) plane of Savage (1993); the function returns where they
%              sit (X,Y).
%
%              Excavation-induced means: start from the laterally constrained
%              gravity state (Sxx = nu/(1-nu)*rg*y, Syy = rg*y) of a flat
%              surface at y = 0, remove the block x < 0, -b < y < 0, and
%              let the new face and lower ground become traction free. The
%              displacement is that of the stress change, i.e. of the
%              stress returned by Savage1993_GravityCliffStress minus the
%              initial state. Plane strain.
%
%              Geometry: top of the cliff y = 0 for x >= 0, cliff face x = 0
%              for -b <= y <= 0, lower ground y = -b for x <= 0.
%
%              REFERENCE: the removed weight is unbounded (the step goes on
%              for ever), so in 2D the displacements grow without bound
%              (like r*log r) and there is no finite 'far' reference.
%              Displacements are given relative to the point (xref, 0) on
%              the top surface: translation and rotation are zero there.
%              Changing xref only adds a rigid-body translation + rotation
%              (strains are unchanged), but that rigid motion grows roughly
%              in proportion to xref: e.g. the cliff-top uplift is about
%              2.2, 24, 240 rg*b^2/mu for xref = 10b, 100b, 1000b (nu = 0.25).
%
%              Savage, W.Z., 1993. Gravity-induced stresses near a vertical
%              cliff. Int. J. Rock Mech. Min. Sci. & Geomech. Abstr. 30(4),
%              325-330.
%
% usage #1:
% [Ux,Uy,X,Y] = Savage1993_GravityCliffDisp(rg,nu,mu,b,u,v)
% [Ux,Uy,X,Y] = Savage1993_GravityCliffDisp(rg,nu,mu,b,u,v,xref)
%
% Arguments: (input)
%       rg    - Weight, the density of material times the acceleration due to
%               gravity.
%
%       nu    - Poisson's Ratio.
%
%       mu    - Shear Modulus.
%
%       b     - Bench (cliff) height.
%
%       u,v   - Observation points in the mapped plane, v <= 0 (see
%               Savage1993_GravityCliffStress).
%
%       xref  - (optional) x of the reference point on the top surface,
%               default 10*b.
%
% Arguments: (output)
%       Ux,Uy - Excavation-induced displacements (x right, y up), relative
%               to the reference point.
%
%       X,Y   - Locations of the observation points.
%
% Method: 2*mu*(Ux + i*Uy) = kappa*phi(Z) - Z*conj(phi'(Z)) - conj(psi(Z)),
%         kappa = 3 - 4*nu, Z = z + i*b, with phi' and psi' the
%         Muskhelishvili potentials of Savage1993_GravityCliffStress (the
%         fitted correction plus the half-plane term H, whose potentials are
%         phi' = i/(2 pi) log(Z), psi' = -i/(2 pi)). phi and psi are
%         integrated numerically (Gauss-Legendre) along straight paths in
%         the w-plane.
%
% Example usage:
%
%  nu = 0.25;
%  mu = 10e9;
%  rg = 2700*9.81;
%  b  = 100;
%  [u,v] = meshgrid(linspace(-6,6,120),-logspace(-4,log10(8),80));
%
%  [Ux,Uy,X,Y] =...
%  Savage1993_GravityCliffDisp(rg,nu,mu,b,u,v);
%
%  quiver(X,Y,Ux,Uy)
%
%  Written with Claude (Anthropic), Nov 2026, for Tim Davis' ElasticSolutions.

if nargin < 7 || isempty(xref), xref = 10*b; end
sz = size(u);
w  = complex(u(:), v(:));
w  = complex(real(w), min(imag(w), -1e-15));
cor = abs(w + 1) < 1e-12;  w(cor) = -1 - 1e-12i;   % base corner: finite displacement
m     = nu/(1-nu);
kappa = 3 - 4*nu;

pk = [ ...
    -1.0000000000000000e+00+3.6631277777468357e-02i
    -1.0000000000000000e+00+1.9205018326236900e-01i
    -1.0000000000000000e+00+6.8478552689813732e-01i
    -1.0000000000000000e+00+2.0000000000000000e+00i
    +1.0000000000000000e+00+3.6631277777468357e-02i
    +1.0000000000000000e+00+1.9205018326236900e-01i
    +1.0000000000000000e+00+6.8478552689813732e-01i
    +1.0000000000000000e+00+2.0000000000000000e+00i
    ];
zeta = [ ...
    -2.2613872124741691e-01+1.5811388300841894e-01i
    +2.5480762113533157e+01+1.4999999999999996e+01i
    -5.0000000000000000e-01+3.1622776601683794e-01i
    -4.9999999999999817e-01+2.9999999999999996e+01i
    -7.7386127875258315e-01+1.5811388300841894e-01i
    -2.6480762113533157e+01+1.4999999999999996e+01i
    ];
lam = [+5.4448373678249995e-01 +9.0852918984599995e-01];
cA = [-7.1765678049440007e-01+2.6829563274373885e-01i, -6.6230562654703423e-02-9.2031814861649064e-01i];
cB = [-2.2503870331787929e-01-6.0194998203343819e-01i, -3.8452860507557540e-01+2.7672545531422031e-02i];
alpha0 = [ ...
    -6.4810986162860382e-05-6.8133162415160795e-05i
    -1.5461080267628730e-03-2.8147939680705412e-03i
    -2.3497297479560589e-02-3.8627950679224185e-02i
    -3.3769797826489262e-01-2.9620810056051416e-01i
    +1.3758529081255821e-05+1.4878761766030264e-06i
    -2.7603842748384150e-04-5.1831188999638196e-04i
    -6.5364685618653712e-03-7.2382302088697678e-03i
    -2.0958959641169509e-02+5.3032685788933047e-03i
    +1.3593739553344007e-02-1.3541910678190749e-02i
    -4.2807773942291197e-04-2.1666798799633213e-04i
    -7.6200807483873922e-02+1.1266254976851338e+00i
    +4.0834508021307513e-03-1.3597456241272465e-03i
    +4.9624022753225319e+00-1.0457000320411140e+01i
    +8.1393259089060168e-05+8.0880822637684744e-04i
    -5.2767699365087370e+00+9.8752046135793741e+00i
    -3.7035058245362067e+00-2.3940330780961352e+00i
    +1.8472861723432263e+01+2.3783895497289041e+01i
    -1.7545575938124955e+01-6.2451564592410733e+01i
    -2.3118397518393952e+00+4.8850375724843964e+00i
    +5.3173696892235816e+00+2.2794395396995770e+00i
    -1.6169739262154014e+01-1.8862533898253638e+01i
    +1.4360399874389543e+01+6.2013069330374442e+01i
    +2.8638155113825943e+00-9.6893020760499198e+00i
    +3.3190652682071486e-01-7.1859065511962539e-01i
    +2.3201796831104088e-01+8.9665853038914645e-01i
    ];
alpha1 = [ ...
    +1.1398301668851489e-04+1.5801647402006194e-04i
    +5.8929159738644148e-03+3.9770896038188062e-03i
    +6.9335103142713086e-02+4.3190433243604463e-02i
    +4.6959577018510357e-01+6.2329874280657817e-02i
    -9.9835920770088358e-06+2.6650083973219736e-05i
    +1.5094381777180238e-03-5.8962208720485200e-04i
    +2.2626949325995450e-02-1.1038150169614581e-02i
    +2.9735366762288610e-02-6.2137436726961037e-02i
    +2.9431248124353916e-02+2.2234465910334927e-02i
    -2.0152383788825406e-05-1.2395674913135504e-03i
    -9.6974088441033879e-01-1.2512795230976561e+00i
    +1.5513206089130268e-04+5.3796718607263764e-03i
    +4.5757965420199245e+00+1.4709380304153163e+01i
    -3.1192133339976271e-04-2.6259424592786976e-04i
    -2.3975664038437050e+00-1.5905569450998728e+01i
    +8.7552102933498555e-01-1.4469871559492982e+00i
    +7.3930505471781620e+00-9.0940393754305386e+00i
    +1.8118946874795583e+01+2.9543418353241815e+01i
    -1.7581084936952773e+00-5.0476785421609787e+00i
    -4.7434114803858496e-01+4.0301032563565427e+00i
    -1.0969140510068048e+01+6.8197936106297092e+00i
    -2.0597758946736235e+01-3.0617156885093749e+01i
    +4.1226764476887050e+00+8.5524353685900003e+00i
    -1.2289824893283459e+00+1.7178724580717103e-01i
    +1.0503892436937534e+00-1.4076669427329638e+00i
    ];
beta0 = [ ...
    -4.2779196020744205e-04-5.8909310826460973e-05i
    -1.2946905099281472e-02+2.2558364182170289e-03i
    -1.3628683261617597e-01+1.6799187107413816e-03i
    -5.6693182393590105e-03-1.4516981927504388e-01i
    +2.9468457755082962e-04-1.7480038797706960e-04i
    -7.8176726609559336e-03-5.5870351071966974e-03i
    -3.6411983108243770e-02-3.3126360188333481e-02i
    -2.9875181616633290e-02-5.6139240637674143e-02i
    +2.0882137827767960e-01+2.7406579172420273e-02i
    +1.5924270627426279e-03-4.5352167947653230e-04i
    -2.1187851379729024e+00+6.1563979676370719e+00i
    -6.2353808145474456e-03+2.3449710252166814e-04i
    +3.8409068662699433e+01-6.0176100168672846e+01i
    +2.7321755446774068e-03-1.8467982670224833e-03i
    -4.5108621997640704e+01+4.7418214611615539e+01i
    +1.7645973747291453e+01-9.5312360631414583e-01i
    -8.2301382252727549e+01+5.9012073257826017e+01i
    +3.7387361314693322e+02+4.4382525681139008e+01i
    -4.6335608662939343e+01-1.2787922884111053e+01i
    -2.7497876990238431e+01+9.4367046069487402e-03i
    +7.5115135112910579e+01-5.3049259074353891e+01i
    -3.8546474448441927e+02-5.6333903889806777e+01i
    +8.1981378093866766e+01+1.9636433190188125e+01i
    +2.0728932979403729e+00-4.2585710520549229e+00i
    +2.1840492840827088e+00+6.3302375126249260e+00i
    ];
beta1 = [ ...
    +9.0369852819936183e-04-5.2389253071635008e-04i
    +1.7452481907240358e-02-2.0985303644215202e-02i
    +1.8082603366415220e-01-1.6145577494113561e-01i
    +1.0916459103715094e-01+5.6674150732768236e-01i
    +1.3809679995348013e-04+3.8761137057690971e-04i
    +1.7950639622485616e-02-1.6883941050504340e-02i
    +1.3174022012022671e-01-7.6868208873750843e-02i
    +2.2977592606007544e-01-8.4581371930812332e-02i
    -6.1860887085533098e-02-4.2263779865366319e-01i
    +3.5846421942674067e-03-1.6877982388878175e-04i
    +1.1573188230671398e+00-7.1698236133298732e+00i
    -7.9107082815709877e-03+6.1620633394359877e-03i
    -3.6007767906683917e+01+1.1118764127137209e+02i
    +1.4447278589825639e-03+3.5715047800007928e-04i
    +5.2622941637251301e+01-9.7207010611017068e+01i
    -2.0888240333089836e+01+2.8324941083920226e+01i
    +8.1362486732611501e+01-2.2062512748689144e+02i
    -7.2752795161411257e+02+4.0947789425117992e+02i
    +9.9220895421333580e+01-2.7279794086982367e+01i
    +3.5200865017637923e+01-4.1419526263960009e+01i
    -6.5789430504050600e+01+1.9322567685209916e+02i
    +7.6338342408081485e+02-3.9815483173602894e+02i
    -1.7321285239701942e+02+5.6991353922331314e+01i
    +1.7853525418615723e+00+8.3001399052801350e+00i
    -1.0085215878989340e+01-6.5179339738668540e+00i
    ];
kappa0 = [+1.0400755474551956e-01 +1.2363228564292161e+00];
kappa1 = [+1.4446485739609832e-01 -1.0933554998992010e+00];

P.alpha = alpha0 + m*alpha1;  P.beta = beta0 + m*beta1;  P.kappa = kappa0 + m*kappa1;
P.pk = pk;  P.zeta = zeta;  P.lam = lam;  P.cA = cA;  P.cB = cB;

% reference point on the top surface: real u > 1 with f(u) = xref/b
wr = refpoint(xref/b);
wa = -2i;                                          % anchor of all integration paths
[phr, psr] = pathint(wa, wr, P);
[Phr, ~, Zr] = pots(wr, P);

% displacement at the points
[ph, ps] = pathint(wa, w, P);
[Ph, ~, Z] = pots(w, P);
U  = kappa*(ph - phr) - Z.*conj(Ph) + Zr*conj(Phr) - conj(ps - psr);   % = 2 mu (Ux + i Uy) / (rg b^2)
U  = U/2;                                          % -> mu (Ux + i Uy)/(rg b^2)
% remove the rigid rotation at the reference point
om = (kappa + 1)*imag(Phr)/2;                      % mu*rotation/(rg b)
U  = U - 1i*om*(Z - Zr);

U  = U*rg*b^2/mu;
Ux = reshape(real(U), sz);  Uy = reshape(imag(U), sz);
X  = reshape(real(Z), sz)*b;  Y = reshape(imag(Z) - 1, sz)*b;
end

% =========================================================================
function [Phi, Psi, Zp, g] = pots(w, P)
% Potentials Phi(z) = phi'(z) and Psi(z) = psi'(z) (frame Z = z + i, b = 1)
% of the excavation field: fitted correction C plus half-plane term H.
w  = w(:);
r2 = sqrt(2);
sm = sqrt(w - 1);  sp = sqrt(w + 1);
Zp = (sm.*sp + 2*(log1pc((w + 1)./((sm - 1i*r2).*(-1i*r2))) + log1pc(sp./sm)))/pi;
z  = Zp - 1i;
g  = (sp./sm)/pi;
za = -1;  zbA = -1i - za;  q = z - za;
Lz = cutlog(Zp) - cutlog(q);
F  = [1./(w - P.pk.'), 1./(z - P.zeta.')];
for e = [0 0.5 1 1.5 2, P.lam, P.lam + 1]
    F = [F, exp(e*Lz)./q];
end
F  = [F, Lz, Lz./q];
Phi = F*P.alpha;  Psi = F*P.beta;
for j = 1:2
    E = exp((P.lam(j) - 1)*(Lz + cutlog(zbA)))*zbA./q;
    Phi = Phi + P.kappa(j)*P.cA(j)*E;  Psi = Psi + P.kappa(j)*P.cB(j)*E;
end
Phi = Phi + 1i/(2*pi)*cutlog(Zp);                  % half-plane term H
Psi = Psi - 1i/(2*pi);
end

function [ph, ps] = pathint(wa, w, P)
% phi(w) - phi(wa) and psi(w) - psi(wa): integrals of Phi dz = Phi f'(w) dw
% along the straight segment wa -> w, split into panels that grow
% geometrically away from wa, the last one graded towards w (corner
% singularities are integrable there).
w = w(:);  n = numel(w);
[xg, wg] = gauss16();
D  = max(abs(w - wa), 1e-300);
K  = max(1, ceil(log2(max(D)/0.5)) + 1);
tb = (0.5./D)*2.^(0:K-1);                          % geometric breaks (n x K)
tb(tb >= 0.999) = 0;
tl = max(tb, [], 2);                               % last break before the end point
tb = [zeros(n,1), tb, tl];
tb = cummax_rows(tb);                              % monotone; unused panels have zero length
ph = zeros(n,1);  ps = ph;
for k = 1:size(tb,2)                               % k = size(tb,2): final panel [tl, 1]
    if k < size(tb,2), t0 = tb(:,k); t1 = tb(:,k+1); last = false;
    else,              t0 = tl;      t1 = ones(n,1); last = true; end
    for j = 1:16
        s = (xg(j) + 1)/2;  ws = wg(j)/2;
        if last                                    % grade towards the end point
            tt = t0 + (t1 - t0).*(1 - (1 - s)^3);  jac = (t1 - t0)*3*(1 - s)^2;
        else
            tt = t0 + (t1 - t0)*s;                 jac = (t1 - t0);
        end
        wq = wa + tt.*(w - wa);
        wq = complex(real(wq), min(imag(wq), -1e-15));
        [Phi, Psi, ~, g] = pots(wq, P);
        dw = (w - wa).*jac*ws;
        ph = ph + Phi.*g.*dw;  ps = ps + Psi.*g.*dw;
    end
end
end

function A = cummax_rows(A)
for k = 2:size(A,2), A(:,k) = max(A(:,k), A(:,k-1)); end
end

function wr = refpoint(xr)
% real u > 1 with (1/pi)[sqrt(u^2-1) + acosh(u)] = xr
wr = max(1 + 1e-9, pi*xr);
for it = 1:100
    f  = (sqrt(wr^2 - 1) + acosh(wr))/pi - xr;
    fp = sqrt((wr + 1)/(wr - 1))/pi;
    wn = max(1 + (wr - 1)/10, wr - f/fp);
    if abs(wn - wr) < 1e-14*wr, wr = wn; break; end
    wr = wn;
end
wr = complex(wr, -1e-15);
end

function [x, w] = gauss16()
x = [-0.9894009349916499; -0.9445750230732326; -0.8656312023878318; -0.7554044083550030; ...
     -0.6178762444026438; -0.4580167776572274; -0.2816035507792589; -0.0950125098376374; ...
      0.0950125098376374;  0.2816035507792589;  0.4580167776572274;  0.6178762444026438; ...
      0.7554044083550030;  0.8656312023878318;  0.9445750230732326;  0.9894009349916499];
w = [0.0271524594117541; 0.0622535239386479; 0.0951585116824928; 0.1246289712555339; ...
     0.1495959888165767; 0.1691565193950025; 0.1826034150449236; 0.1894506104550685; ...
     0.1894506104550685; 0.1826034150449236; 0.1691565193950025; 0.1495959888165767; ...
     0.1246289712555339; 0.0951585116824928; 0.0622535239386479; 0.0271524594117541];
end

function L = cutlog(q)
a = angle(q);  a(a > 0.75*pi) = a(a > 0.75*pi) - 2*pi;
L = log(abs(q)) + 1i*a;
end

function y = log1pc(x)
u = 1 + x;  y = log(u);  k = (u ~= 1);
y(k) = y(k).*x(k)./(u(k) - 1);  y(~k) = x(~k);
end
