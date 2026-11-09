% Savage1993_GravityCliffStress_Plots
% Script: gravity stresses beneath a vertical cliff (bench height b), drawn
% like Figs 1, 3, 4 and 5 of Savage (1993). It calls both functions:
%   Savage1993_GravityCliffStress         (mapped u,v plane -> curvilinear grid)
%   Savage1993_GravityCliffStress_CoordsIn (real x,y grid)
% and overlays the two as a check (black = mapped grid, red dashed = x,y grid).
% A further figure shows the stresses along the free surface. If doDisp is
% true the excavation-induced displacements are drawn too, using
%   Savage1993_GravityCliffDisp and Savage1993_GravityCliffDisp_CoordsIn
% (relative to the point (xref, 0) on the top surface).
% Stresses are shown as sigma/(rg*b), tension positive; axes are x/b, y/b.
%
% Savage, W.Z., 1993. Gravity-induced stresses near a vertical cliff.
% Int. J. Rock Mech. Min. Sci. & Geomech. Abstr. 30(4), 325-330.

clear; close all

% ---------------- inputs ----------------
rg = 2700*9.81;      % weight, rho*g
nu = 1/3;            % Poisson's ratio (Savage's figures use 1/3)
mu = 10e9;           % shear modulus (not used by the stresses)
b  = 100;            % bench (cliff) height
xlim_b = [-1 1.5];   % plotting window in units of b (as in the paper)
ylim_b = [-3 0];
saveFigs = false;    % true -> write Savage1993_Fig*.png
doDisp   = true;     % also draw the excavation-induced displacements
xref     = 10*b;     % displacement reference point (xref, 0) on the top surface
exag     = 0.3;      % plot exaggeration: largest displacement drawn as exag*b

% ---------------- 1) mapped (u,v) plane ----------------
% grid clustered at the corners (u = -1, 1) and at the surface (v = 0)
du = logspace(-5, log10(11), 150);
u  = unique([-1 - du, -1 + du(du < 1), 1 - du(du < 1), 1 + du, linspace(-1, 1, 101)]);
v  = -logspace(-5, log10(12), 300);
[U, V] = meshgrid(u, v);
[Sxx, Syy, Sxy, X, Y] = Savage1993_GravityCliffStress(rg, nu, mu, b, U, V);
c  = (Sxx + Syy)/2;  r = sqrt(((Sxx - Syy)/2).^2 + Sxy.^2);
Fmap = {(c + r)/(rg*b), (c - r)/(rg*b), r/(rg*b)};          % s1, s3, tmax

% ---------------- 2) real (x,y) coordinates ----------------
[xg, yg] = meshgrid(linspace(xlim_b(1), xlim_b(2), 126)*b, linspace(ylim_b(1), ylim_b(2), 151)*b);
[Sxx2, Syy2, Sxy2] = Savage1993_GravityCliffStress_CoordsIn(rg, nu, mu, b, xg, yg);
c2 = (Sxx2 + Syy2)/2;  r2 = sqrt(((Sxx2 - Syy2)/2).^2 + Sxy2.^2);
Fxy = {(c2 + r2)/(rg*b), (c2 - r2)/(rg*b), r2/(rg*b)};

% ---------------- Fig 1: the curvilinear coordinates ----------------
figure('Color', 'w', 'Name', 'Fig 1'); hold on
ul = [-40 -10 -5 -3 -2 -1.5 -1.2 -1 -0.8 -0.5 -0.25 0 0.25 0.5 0.8 1 1.2 1.5 2 3 5 10 40];
vl = -[0.05 0.1 0.2 0.35 0.5 0.75 1 1.5 2 3 5 10];
ut = [linspace(-60, -1.0001, 400), linspace(-1, 1, 200), linspace(1.0001, 60, 400)];
vt = -[logspace(-4, 0, 150), linspace(1.02, 40, 200)];
for k = 1:numel(vl)
    [~, ~, ~, Xl, Yl] = Savage1993_GravityCliffStress(rg, nu, mu, b, ut, vl(k) + 0*ut);
    plot(Xl/b, Yl/b, 'Color', [0.6 0.6 0.6]);
end
for k = 1:numel(ul)
    [~, ~, ~, Xl, Yl] = Savage1993_GravityCliffStress(rg, nu, mu, b, ul(k) + 0*vt, vt);
    plot(Xl/b, Yl/b, 'Color', [0.6 0.6 0.6]);
end
plot([-2 0 0 2], [-1 -1 0 0], 'k', 'LineWidth', 2);
axis equal; axis([-2 2 -3 0.5]); box on
xlabel('x/b'); ylabel('y/b'); title('Fig 1: lines of constant u and v');

% ---------------- Figs 3-5: s1, s3, tmax ----------------
names = {'\sigma_1/(\rho g b)', '\sigma_3/(\rho g b)', '\tau_{max}/(\rho g b)'};
dls   = [0.2 0.2 0.1];                      % contour intervals as in the paper
for f = 1:3
    figure('Color', 'w', 'Name', sprintf('Fig %d', f + 2), 'Position', [100 100 520 640]); hold on
    dl = dls(f);
    lv = round((-4:dl:4)/dl)*dl;
    heavy = mod(round(lv/dl), 2) == 0;      % every second contour heavy and labelled
    [~, h1] = contour(X/b, Y/b, Fmap{f}, lv(~heavy), 'k');  set(h1, 'LineWidth', 0.5);
    [C2, h2] = contour(X/b, Y/b, Fmap{f}, lv(heavy), 'k');  set(h2, 'LineWidth', 1.6);
    try, clabel(C2, h2, 'FontSize', 8, 'LabelSpacing', 400); catch, end
    [~, h3] = contour(xg/b, yg/b, Fxy{f}, lv, 'r--');        set(h3, 'LineWidth', 0.8);
    fill([xlim_b(1) 0 0 xlim_b(1)], [-1 -1 0 0], 'w', 'EdgeColor', 'none');   % air
    plot([xlim_b(1) 0 0 xlim_b(2)], [-1 -1 0 0], 'k', 'LineWidth', 1.5);
    axis equal; axis([xlim_b ylim_b]); box on
    set(gca, 'XAxisLocation', 'top', 'TickDir', 'in');
    xlabel('x/b'); ylabel('y/b');
    title(sprintf('Fig %d: %s,  \\nu = %.3g (tension positive)', f + 2, names{f}, nu));
    if saveFigs, print('-dpng', '-r150', sprintf('Savage1993_Fig%d.png', f + 2)); end
end
if saveFigs, figure(1); print('-dpng', '-r150', 'Savage1993_Fig1.png'); end

% ---------------- stresses along the free surface ----------------
% Walk along the ground surface (v = 0 in the mapped plane, u from +inf to
% -inf): top of the cliff -> top corner -> down the face -> base corner ->
% lower ground. s = distance along the free surface from the cliff top
% corner: positive along the top surface, 0..-b down the face, < -b along
% the lower ground (so the cliff top is on the right, the base on the left).
e  = logspace(-6, 0, 200);
us = unique([1 + e, linspace(2, 12, 100), 1 - e, -1 + e, linspace(-0.99, 0.99, 100), -1 - e, -linspace(2, 30, 150)]);
us = fliplr(us);                                 % u decreasing = moving along the surface
vs = zeros(size(us));                            % on the surface (limit from below)
[Sxxs, Syys, Sxys, Xs, Ys] = Savage1993_GravityCliffStress(rg, nu, mu, b, us, vs);
sdist = zeros(size(us));
top  = us > 1;  face = abs(us) <= 1;  low = us < -1;
sdist(top)  =  Xs(top);                          % along the top surface
sdist(face) =  Ys(face);                         % down the face
sdist(low)  =  Xs(low) - b;                      % along the lower ground
cs = (Sxxs + Syys)/2;  rs = sqrt(((Sxxs - Syys)/2).^2 + Sxys.^2);
S1s = cs + rs;  S2s = cs - rs;
keep = sdist > -4*b & sdist < 2*b;
sb = sdist(keep)/b;  n1 = rg*b;
figure('Color', 'w', 'Name', 'Free surface', 'Position', [100 100 900 650]);
lims = [-4 2];
away = abs(sb + 1) > 0.05;                       % y-limits ignore the singular base corner
for p = 1:2
    subplot(2, 1, p); hold on; box on
    if p == 1
        plot(sb, Sxxs(keep)/n1, 'b', 'LineWidth', 1.4);
        plot(sb, Syys(keep)/n1, 'r', 'LineWidth', 1.4);
        plot(sb, Sxys(keep)/n1, 'k', 'LineWidth', 1.4);
        F = [Sxxs(keep); Syys(keep); Sxys(keep)]/n1;
        legend({'\sigma_{xx}', '\sigma_{yy}', '\sigma_{xy}'}, 'Location', 'eastoutside');
        title(sprintf('Stresses along the free surface, \\nu = %.3g (tension positive)', nu));
    else
        plot(sb, S1s(keep)/n1, 'b', 'LineWidth', 1.4);
        plot(sb, S2s(keep)/n1, 'r', 'LineWidth', 1.4);
        F = [S1s(keep); S2s(keep)]/n1;
        legend({'\sigma_1', '\sigma_2'}, 'Location', 'eastoutside');
    end
    F = F(:, away);  yl = [min(F(:)) max(F(:))];  yl = yl + [-0.1 0.1]*diff(yl);
    plot(lims, [0 0], 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    plot([0 0], yl, 'k--', 'HandleVisibility', 'off');  plot([-1 -1], yl, 'k--', 'HandleVisibility', 'off');
    text(0, yl(2), ' cliff top', 'VerticalAlignment', 'top');
    text(-1, yl(2), ' base corner', 'VerticalAlignment', 'top');
    text(1, yl(1), 'top surface', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(-0.5, yl(1), 'face', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(-2.5, yl(1), 'lower ground', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    xlim(lims); ylim(yl);
    ylabel('\sigma/(\rho g b)');
    if p == 2, xlabel('s/b, distance along the free-surface from the cliff top (negative towards the base)'); end
end
if saveFigs, print('-dpng', '-r150', 'Savage1993_FreeSurface.png'); end

% ---------------- excavation-induced displacements ----------------
if doDisp
    nd = rg*b^2/mu;                                % displacement scale: plotted u/nd, real u = (plotted)*rg*b^2/mu
    % (a) field: |u| contours on the mapped grid, arrows from the x,y grid
    Ud = U(1:2:end, 1:3:end);  Vd = V(1:2:end, 1:3:end);     % coarser grid is plenty here
    [Uxm, Uym, Xm, Ym] = Savage1993_GravityCliffDisp(rg, nu, mu, b, Ud, Vd, xref);
    [xq, yq] = meshgrid(linspace(xlim_b(1), xlim_b(2), 14)*b, linspace(ylim_b(1), ylim_b(2), 16)*b);
    [Uxq, Uyq] = Savage1993_GravityCliffDisp_CoordsIn(rg, nu, mu, b, xq, yq, xref);
    figure('Color', 'w', 'Name', 'Displacement', 'Position', [100 100 560 640]); hold on
    Um = hypot(Uxm, Uym)/nd;
    inwin = Xm/b >= xlim_b(1) & Xm/b <= xlim_b(2) & Ym/b >= ylim_b(1);
    lv = linspace(min(Um(inwin)), max(Um(inwin)), 15);
    [Cd, hd] = contour(Xm/b, Ym/b, Um, lv);
    try, set(hd, 'LineColor', [0.5 0.5 0.5]); catch, end
    try, clabel(Cd, hd, 'FontSize', 7, 'LabelSpacing', 500); catch, end
    quiver(xq/b, yq/b, Uxq/nd, Uyq/nd, 1.5, 'k');
    % deformed surface, exaggerated
    us_ = fliplr(unique([1 + logspace(-6, 0, 100), linspace(2, 12, 60), 1 - logspace(-6, 0, 60), ...
                -1 + logspace(-6, 0, 60), linspace(-0.99, 0.99, 60), -1 - logspace(-6, 0, 100), -linspace(2, 30, 80)]));
    [Uxs, Uys, Xs2, Ys2] = Savage1993_GravityCliffDisp(rg, nu, mu, b, us_, 0*us_, xref);
    sc = exag*b/max(hypot(Uxs, Uys));
    fill([xlim_b(1) 0 0 xlim_b(1)], [-1 -1 0 0], 'w', 'EdgeColor', 'none');
    plot([xlim_b(1) 0 0 xlim_b(2)], [-1 -1 0 0], 'k', 'LineWidth', 1.5);
    plot((Xs2 + sc*Uxs)/b, (Ys2 + sc*Uys)/b, 'r', 'LineWidth', 1.5);
    plot(xref/b, 0, 'kv', 'MarkerFaceColor', 'k');
    axis equal; axis([xlim_b ylim_b(1) 0.6]); box on
    set(gca, 'XAxisLocation', 'top', 'TickDir', 'in');
    xlabel('x/b'); ylabel('y/b');
    title({'Contours and arrows: |u|/(\rho g b^2/\mu)', sprintf('Red: deformed surface, u exaggerated \\times%.3g;  reference at x/b = %.3g', sc, xref/b)});
    if saveFigs, print('-dpng', '-r150', 'Savage1993_Displacement.png'); end

    % (b) along the free surface
    [Uxs, Uys, Xs3, Ys3] = Savage1993_GravityCliffDisp(rg, nu, mu, b, us, vs, xref);
    keep = sdist > -4*b & sdist < 2*b;  sb = sdist(keep)/b;
    figure('Color', 'w', 'Name', 'Free surface displacement', 'Position', [100 100 900 420]); hold on; box on
    plot(sb, Uxs(keep)/nd, 'b', 'LineWidth', 1.4);
    plot(sb, Uys(keep)/nd, 'r', 'LineWidth', 1.4);
    legend({'u_x', 'u_y'}, 'Location', 'eastoutside');
    F = [Uxs(keep); Uys(keep)]/nd;  yl = [min(F(:)) max(F(:))];  yl = yl + [-0.1 0.1]*diff(yl);
    plot([-4 2], [0 0], 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    plot([0 0], yl, 'k--', 'HandleVisibility', 'off');  plot([-1 -1], yl, 'k--', 'HandleVisibility', 'off');
    text(0, yl(2), ' cliff top', 'VerticalAlignment', 'top');
    text(-1, yl(2), ' base corner', 'VerticalAlignment', 'top');
    text(1, yl(1), 'top surface', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(-0.5, yl(1), 'face', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    text(-2.5, yl(1), 'lower ground', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    xlim([-4 2]); ylim(yl);
    xlabel('s/b, distance along the free-surface from the cliff top (negative towards the base)');
    ylabel('u/(\rho g b^2/\mu)');
    title(sprintf('Excavation displacement along the free surface (reference at x/b = %.3g), \\nu = %.3g', xref/b, nu));
    if saveFigs, print('-dpng', '-r150', 'Savage1993_FreeSurfaceDisplacement.png'); end
end


% ---------------- check: the two functions agree ----------------
[~, ~, ~, Xc, Yc] = Savage1993_GravityCliffStress(rg, nu, mu, b, U(1:7:end, 1:7:end), V(1:7:end, 1:7:end));
[a1, a2, a3] = Savage1993_GravityCliffStress(rg, nu, mu, b, U(1:7:end, 1:7:end), V(1:7:end, 1:7:end));
[b1, b2, b3] = Savage1993_GravityCliffStress_CoordsIn(rg, nu, mu, b, Xc, Yc);
d = max(abs([a1(:) - b1(:), a2(:) - b2(:), a3(:) - b3(:)]), [], 2);
s = max(rg*b, max(abs([a1(:), a2(:), a3(:)]), [], 2));     % relative near the singular corner
fprintf('max relative difference, mapped vs CoordsIn = %.1e\n', max(d(isfinite(d))./s(isfinite(d))));