function [CriticalVolume] = CriticalVolumeDavis2020(nu,mu,Kc,delta_gamma)
% CriticalVolumeDavis2020: A volume of a penny-shaped fracture such that
% the upper tip is K+=Kc and the lower K-=0. Eq.6 Davis et al 2020:
% Critical Fluid Injection Volumes for Uncontrolled Fracture Ascent
% Geophysical research letters
%
% Arguments: (input)
%     nu            - Poisson's ratio of host material.
%
%     mu            - Shear modulus of the host material.
%
%     Kc            - Fracture toughness of the host material.
%
%     delta_gamma   - Stress gradient acting on the fractures walls, for a
%                     vertical fracture: delta_gamma=(rho_r-rho_f)g. Where
%                     rho_r is the rock density, rho_f the fluid density
%                     and g the gradient in stress.
%
% Arguments: (output)
% 	  CriticalVolume  - The theoretical volume of the fracture before this
%                       begins to ascend upwards through the crust.%
%
%  Author: Tim Davis

CriticalVolume=(((1-nu)/(16*mu))*((9*pi^4*Kc^8)/(delta_gamma^5))^(1/3));

end

