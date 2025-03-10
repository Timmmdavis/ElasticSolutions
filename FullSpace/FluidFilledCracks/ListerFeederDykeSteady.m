% Script to calculate dyke width w(x,z) and lateral extent b(z) 
% Based on equations 20a and 20b from Lister (1991)
% "Steady solutions for feeder dykes in a density-stratified lithosphere"
% Lateral extent limited by viscous flow (Sec.3.4)

clear
close all

% Physical constants
g = 9.81;            % Gravitational acceleration (m/s^2)
eta = 100;           % Magma viscosity (Pa s)
m = 20e9;            % Elastic modulus = G/(1-v) where G is shear modulus (Pa)
q = 100;             % Volumetric flow rate (m^3/s)

% Density parameters
rho_s0 = 3000;       % Reference solid density at z=0 (kg/m^3) - goes to fluid at level of neutral bouy
rho_f = 2700;        % Fluid (magma) density (kg/m^3)
z_L = 20000;         % Level of neutral buoyancy (m)

%Calculations
drho_0 = rho_s0-rho_f;                     % Initial density difference at z=0 (kg/m^3)
R=drho_0/z_L; %Change in rock density with depth (kg/m^2)

% Calculate density difference profile
z_fine = linspace(0, z_L, 5000);  % Fine depth array for integration (m)
z = linspace(0, z_L, 500);        % Depth array for plotting (m)
drho_fine = drho_0 * (1 - z_fine/z_L);  % Linear decrease in density difference (fine grid)
drho = drho_0 * (1 - z/z_L);            % Linear decrease in density difference (plotting grid)

% Calculate z_hat (equation 18)
% z_hat(z) = ∫ [Δρ(0)/Δρ(ξ)]^(4/3) dξ
z_hat = zeros(size(z));
for i = 1:length(z)
    % Integration from 0 to z(i)
    idx = find(z_fine <= z(i));
    if ~isempty(idx)
        integrand = (drho_0./drho_fine(idx)).^(4/3);
        dz = diff(z_fine(idx([1, end]))) / (length(idx) - 1);
        z_hat(i) = sum(integrand) * dz;
    end
end

% Calculate lateral extent b(z) using equation 20b
% b(z) = 2.62((q*η*m^2)/(g*Δρ(0))^4)^(1/10) * ẑ^(1/10)
b = zeros(size(z));
for i = 1:length(z)
    b(i) = 2.62 * ((q*eta*m^3*z_hat(i)^3) / (g*drho_0)^4)^(1/10) ;
end
b(1)=nan;
b(end)=nan;

% Calculate width w(x,z) using equation 20a
% w(x,z) = 0.904(Δρ(0)/Δρ(z))^(1/3) * ((q^3*η^3)/(m*(g*Δρ(0))^2*ẑ^2))^(1/10) * (1-(x/b)^2)^(3/2)
x = linspace(-nanmax(b)*0.95, nanmax(b)*0.95, 50);  % Avoid exact boundaries for numerical stability
[X, Z] = meshgrid(x, z);
W = zeros(size(X));

for i = 1:length(z)
    for j = 1:length(x)
        if abs(x(j)) < b(i)  % Stay within lateral extent
            W(i,j) = 0.904 .* (drho_0/drho(i)).^(1/3) * ...
                    ((q^3*eta^3)/(m*(g*drho_0)^2*z_hat(i)))^(1/10) .* ...
                    (1 - (x(j).^2./b(i).^2)).^(3/2);
        end
    end
end

% Plot results
figure(1)
plot(b/1000, z/1000, 'LineWidth', 2)
%set(gca, 'YDir', 'reverse')  % Invert y-axis to show depth increasing downward
xlabel('Lateral extent b(z) (km)')
ylabel('Depth z (km)')
title('Dyke Lateral Extent vs Depth')
grid on

figure(2);hold on
h=pcolor(X/1000, Z/1000, W*1000) ; % Convert to km and mm
shading flat;%remove edges
contour(X/1000, Z/1000, W*1000,'k','ShowText','on')
%set(gca, 'YDir', 'reverse')  % Invert y-axis to show depth increasing downward
colorbar
xlabel('Horizontal position x (km)')
ylabel('Depth z (km)')
title('Dyke Width w(x,z) (mm)')

% Plot width profile at different depths
figure(3)
depths_to_plot = [0.1, 0.5, 0.9] * z_L;  % Plot at 10%, 50%, and 90% of the LNB depth
colors = {'r-', 'g-', 'b-'};
legend_entries = cell(length(depths_to_plot), 1);

hold on
for i = 1:length(depths_to_plot)
    depth_idx = find(abs(z - depths_to_plot(i)) == min(abs(z - depths_to_plot(i))), 1);
    plot(x/1000, W(depth_idx,:)*1000, colors{i}, 'LineWidth', 2)
    legend_entries{i} = sprintf('z = %.1f km', z(depth_idx)/1000);
end
hold off
xlabel('Horizontal position x (km)')
ylabel('Dyke width w (mm)')
title('Width Profiles at Different Depths')
legend(legend_entries)
grid on

% Print key values
fprintf('Maximum lateral extent: %.2f km\n', max(b)/1000)
fprintf('Width at center (x=0), z=0: %.2f mm\n', W(1, ceil(length(x)/2))*1000)
fprintf('Width at center (x=0), z=%.1f km (near LNB): %.2f mm\n', ...
        z(end)/1000, W(end, ceil(length(x)/2))*1000)
        
% Show analytical expressions
fprintf('\nAnalytical expressions from the paper:\n')
fprintf('Equation 20a: w(x,z) = 0.904(Δρ(0)/Δρ(z))^(1/3) * ((q^3*η^3)/(m*(g*Δρ(0))^2*ẑ^2))^(1/10) * (1-(x/b)^2)^(3/2)\n')
fprintf('Equation 20b: b(z) = 2.62((q*η*m^2)/(g*Δρ(0))^4)^(1/10) * ẑ^(1/10)\n')
fprintf('Where ẑ(z) = ∫[Δρ(0)/Δρ(ξ)]^(4/3) dξ\n')