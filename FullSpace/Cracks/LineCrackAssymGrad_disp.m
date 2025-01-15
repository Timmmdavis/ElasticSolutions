function [Dn,Area] = LineCrackAssymGrad_disp(c,gamma,z,Eprime)
% LineCrackAssymGrad_disp: Wall openings of a line crack line crack
% with an hourglass shaped stress gradient (linear) that reaches gamma*c at
% the upper tip -gamma*c at the lower and gradient is 0 at the crack
% centre.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     Eprime       - Elastic constant of the material surronding the fracture
%                    E/(1-nu^2).
%
%     gamma        - Linear stress gradient acting on the crack walls.
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
% 	    gamma*c
%       ――――
% \←---|---→/
%  \←--|--→/      z
%   \←-|-→/       ↑
%    \←|→/        | 
%     \|/          --→x   
%     /|\     |
%    /→|←\    |
%   /-→|←-\   |c
%  /--→|←--\  |
% /---→|←---\ |
% 	    ――――	
% 	    -gamma*c
%
% Example usage:
% c=1;
% E=1e9; %Youngs Mod
% nu=0.25; %Poisson's ratio
% Eprime=E/(1-nu^2);
% gamma=1000*9.81;
% z=linspace(-c,c,120)';
% x=zeros(size(z));
% [Dn,Area] = LineCrackAssymGradient_disp(c,gamma,z,Eprime);
% zwall_grad=[z;flipud(z)]; %merging pos & neg side
% xwall_grad=[x+Dn/2;flipud(x-Dn/2)];
% figure;
% patch(xwall_grad,zwall_grad,[0.9,0.9,0.9]);
% hold on;%grid on
% scatter(x+Dn/2,z,'r.');  %pos side
% scatter(x-Dn/2,z,'b.');  %neg side
%
%  Author: Tim Davis

%Tada P.154
Dn=((2*(gamma*c)*c)/Eprime).*(z./c).*sqrt(1-(z/c).^2);

Area=0.;

end

