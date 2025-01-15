function [Dn,Area] = LineCrackConstPressure_disp(c,p,Eprime,z)
% LineCrackConstPressure_disp: Wall openings of a line crack 
% with a constant pressure across its faces.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     Eprime       - Elastic constant of the material surronding the fracture
%                    E/(1-nu^2).
%
%     p            - Constant pressure
%
%     z            - Location on wall above crack centre.
%
% Arguments: (output)
% 	  Dn           - The seperation of the fracture's walls.
%
%     Area         - The area of the gap opened between the fracture's
%                    walls (whole fracture).
%
% Diagram:
%
% 	     p
%       ――――
% |←---|---→|
% |←---|---→|       z
% |←---|---→|       ↑
% |←---|---→|       |
% |←---|---→|        --→x  
% |←---|---→| |
% |←---|---→| |
% |←---|---→| |c
% |←---|---→| |
% |←---|---→| |
%
% Example usage:
% c=1;
% E=1e9; %Youngs Mod
% nu=0.25; %Poisson's ratio
% Eprime=E/(1-nu^2);
% p=1e6;
% z=linspace(-c,c,120)';
% x=zeros(size(z));
% [Dn,Area] = LineCrackConstPressure_disp(c,p,Eprime,z);
% zwall_grad=[z;flipud(z)]; %merging pos & neg side
% xwall_grad=[x+Dn/2;flipud(x-Dn/2)];
% figure;
% patch(xwall_grad,zwall_grad,[0.9,0.9,0.9]);
% hold on;%grid on
% scatter(x+Dn/2,z,'r.');  %pos side
% scatter(x-Dn/2,z,'b.');  %neg side
%
%  Author: Tim Davis

%Tada P.125
Dn=((4*p)/(Eprime))*sqrt(c^2-z.^2);

Area=(2*p*pi*c^2)/Eprime;

end

