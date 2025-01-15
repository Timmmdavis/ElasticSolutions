% CompareAnalyticalAndNumericalCrackSolutions.m
% Script to compare different analytical and numerical solutions for crack problems
% Based on Tada's Handbook and Townsend & Pollard (2018)

% Clear workspace and figures
clear
close all

% Define basic parameters
c = 15;              % Crack half-length
n = 150;             % Number of discretization points
Eprime = 1e9;       % Plane strain elastic modulus
nu = 0.1;          % Poisson's ratio
mu = Eprime/(2*(1+nu)); % Shear modulus
Rg = 1e3;          % Gravitational term
Kc = 1e5;          % Fracture toughness
pf = (Rg*(c/2)^2)+ Kc/sqrt(pi*c);   % Fluid pressure
wf = 0.001;        % Reference crack width
IntRelTol = 1e-6;  % Integration relative tolerance
IntAbsTol = 1e-9;  % Integration absolute tolerance

% Create position arrays
z = linspace(-c, c, n);
zmid = (z(1:end-1) + z(2:end))/2;

% Precompute parameters for Pollard-Townsend methods
Precomp = TownsendPollardPrecompute2018(c, n, nu);

% Initialize figure
figure('Position', [100 100 1200 1000])

% Case A: Quadratic pressure distribution
subplot(3,2,1)
a_TnFunc = @(z) pf - Rg.*z.^2/2;
a_wfan = ((1/6.*Rg.*(c^2-z.^2)+Kc./sqrt(c*pi)).*(4./Eprime.*sqrt(c^2-z.^2)));
a_Kcan = Kc;
a_Aan = (8*Kc*c^(3/2)*pi^(1/2)+Rg*c^4*pi)/(4*Eprime);
[a_KuPT, a_KlPT, a_wfPT, a_APT] = PollardTractionIntegrationVerticalCrack(a_TnFunc, c, mu, nu, z, zmid, Precomp);
[a_KuPTL, a_KlPTL, a_wfPTL, a_APTL] = PollardTractionIntegrationVerticalCrackWithLinGrads(a_TnFunc, c, mu, nu, z, zmid, Precomp);
[a_KuTada, a_KlTada, a_wfTada, a_ATada] = TadaTractionIntegrationVerticalCrack(a_TnFunc, c, z, Eprime, IntRelTol, IntAbsTol);
PlotCrackComparison(z, a_wfan, a_wfPT, a_wfPTL, a_wfTada, 'Quadratic Pressure')

% Case B: Constant pressure
subplot(3,2,2)
b_TnFunc = @(z) pf + 0*z;  % Constant function
b_wfan = ((pf).*(4./Eprime.*sqrt(c^2-z.^2)));
b_Kcan = pf*sqrt(c*pi);
b_Aan = (2*pf*c^(2)*pi)/(Eprime);
[b_KuPT, b_KlPT, b_wfPT, b_APT] = PollardTractionIntegrationVerticalCrack(b_TnFunc, c, mu, nu, z, zmid, Precomp);
[b_KuPTL, b_KlPTL, b_wfPTL, b_APTL] = PollardTractionIntegrationVerticalCrackWithLinGrads(b_TnFunc, c, mu, nu, z, zmid, Precomp);
[b_KuTada, b_KlTada, b_wfTada, b_ATada] = TadaTractionIntegrationVerticalCrack(b_TnFunc, c, z, Eprime, IntRelTol, IntAbsTol);
PlotCrackComparison(z, b_wfan, b_wfPT, b_wfPTL, b_wfTada, 'Constant Pressure')

% Case C: Constant displacement
subplot(3,2,3)
c_TnFunc = @(z) -(wf*Eprime*c)./(2*pi*(z.^2 - c^2));
c_wfan = ones(size(z))*wf;
c_Kcan = nan;
c_Aan = 2*c*wf;
[c_KuPT, c_KlPT, c_wfPT, c_APT] = PollardTractionIntegrationVerticalCrack(c_TnFunc, c, mu, nu, z, zmid, Precomp);
[c_KuPTL, c_KlPTL, c_wfPTL, c_APTL] = PollardTractionIntegrationVerticalCrackWithLinGrads(c_TnFunc, c, mu, nu, z, zmid, Precomp);
[c_KuTada, c_KlTada, c_wfTada, c_ATada] = TadaTractionIntegrationVerticalCrack(c_TnFunc, c, z, Eprime, IntRelTol, IntAbsTol);
PlotCrackComparison(z, c_wfan, c_wfPT, c_wfPTL, c_wfTada, 'Constant Displacement')

% Case D: Linear decrease to tips
subplot(3,2,4)
d_TnFunc = @(z) pf*(1-abs(z)/c);
d_wfan = (4.*pf.*c)./Eprime.*((1-1./pi).*sqrt(1-(z./c).^2)-1./pi.*(z./c).^2.*acosh(c./abs(z)));
d_Kcan = (1-2/pi).*pf.*sqrt(pi.*c);
d_Aan = 2*(pi-4/3)*(pf*c^2)/Eprime;
[d_KuPT, d_KlPT, d_wfPT, d_APT] = PollardTractionIntegrationVerticalCrack(d_TnFunc, c, mu, nu, z, zmid, Precomp);
[d_KuPTL, d_KlPTL, d_wfPTL, d_APTL] = PollardTractionIntegrationVerticalCrackWithLinGrads(d_TnFunc, c, mu, nu, z, zmid, Precomp);
[d_KuTada, d_KlTada, d_wfTada, d_ATada] = TadaTractionIntegrationVerticalCrack(d_TnFunc, c, z, Eprime, IntRelTol, IntAbsTol);
PlotCrackComparison(z, d_wfan, d_wfPT, d_wfPTL, d_wfTada, 'Linear Decrease to Tips')

% Case E: Linear increase to tips
subplot(3,2,5)
e_TnFunc = @(z) pf*(abs(z)/c);
e_wfan = (4.*pf.*c)./(pi*Eprime).*(sqrt(1-(z./c).^2)+(z./c).^2.*acosh(c./abs(z)));
e_Kcan = 2/pi.*pf.*sqrt(pi.*c);
e_Aan = (8/3)*(pf*c^2)/Eprime;
[e_KuPT, e_KlPT, e_wfPT, e_APT] = PollardTractionIntegrationVerticalCrack(e_TnFunc, c, mu, nu, z, zmid, Precomp);
[e_KuPTL, e_KlPTL, e_wfPTL, e_APTL] = PollardTractionIntegrationVerticalCrackWithLinGrads(e_TnFunc, c, mu, nu, z, zmid, Precomp);
[e_KuTada, e_KlTada, e_wfTada, e_ATada] = TadaTractionIntegrationVerticalCrack(e_TnFunc, c, z, Eprime, IntRelTol, IntAbsTol);
PlotCrackComparison(z, e_wfan, e_wfPT, e_wfPTL, e_wfTada, 'Linear Increase to Tips')

% Add overall title
sgtitle('Comparison of Crack Opening Solutions')

% Helper function for consistent plotting
function PlotCrackComparison(z, wfan, wfPT, wfPTL, wfTada, titleStr)
    hold on
    plot(z, wfan, 'k-', 'LineWidth', 2, 'DisplayName', 'Analytical')
    plot(z, wfPT, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Pollard-Townsend')
    plot(z, wfPTL, 'b-.', 'LineWidth', 1.5, 'DisplayName', 'PT with Linear')
    plot(z, wfTada, 'g:', 'LineWidth', 1.5, 'DisplayName', 'Tada Integration')
    grid on
    xlabel('Position along crack (z/c)')
    ylabel('Opening width')
    title(titleStr)
    legend('show', 'Location', 'best')
    hold off
end