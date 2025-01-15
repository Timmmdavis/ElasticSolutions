function [Dn]=PennyCrackAssymGrad_disp(y,z,gamma,mu,nu,c)
% PennyCrackAssymGrad_disp: Compute the separation of the cracks faces
% (D_n) due to linear stress gradient gamma.
%
% Simplified opening from elliptical equations (Tim Davis)
%
% Arguments: (input)
%     c             - Crack radius.
%
%     nu            - Poisson's ratio of host material.
%
%     mu            - Shear modulus of the host material.
%
%     y             - Horizontal coordinate on the crack wall, relative to
%                     the centre (cartesian coords, on the crack).
%
%     z             - Vertical coordinate on the crack wall, relative to
%                     the centre (cartesian coords, on the crack).
%
%     gamma         - Linear stress gradient, the stress gradient acts along
%                     the z-axis, is asymmetric, 0 at the centre, magnitude
%                     gamma*c at the top tip and -gamma*c at the basal
%                     tip.
%
% Arguments: (output)
%     Dn          - Seperation of the faces due to linear stress gradient
%                     gamma.

% Diagram:
% 1. Looking at crack face (yz): 
% 
%                   ooo OOO OOO ooo
%               oOO        |        OOo
%           oOO            |            OOo
%        oOO               |c (z)         OOo
%      oOO                 |                 OOo
%    oOO                   | Phi (Φ)          OOo
%   oOO                    |   /               OOo
%  oOO                     |__/                 OOo
% oOO                      | /                   OOo
% oOO                      |/_____________________OOo
% oOO                             c (y)           OOo
% oOO                                             OOo
% oOO                                             OOo
%  oOO                                           OOo
%   oOO                                         OOo
%    oOO                                       OOo
%      oOO                                   OOo
%        oO                                OOo
%           oOO                         OOo
% z-axis        oOO                 OOo
% ^                 ooo OOO OOO ooo
% |
% |
% |
% |
% --------------> y-axis
%
% 2. Assymetric gradient boundary condition, 2D cross section (xz):
%
% 	   
% 	    gamma*c
%       ――――	
% \←---|---→/
%  \←--|--→/
%   \←-|-→/
%    \←|→/
%     \|/
%     /|\     |
%    /→|←\    |
%   /-→|←-\   |c (z-axis)
%  /--→|←--\  |
% /---→|←---\ |
% 	    ――――	
% 	   -gamma*c
%
%
% Example usage:
% gamma=1000*9.81;
% c=500;
% smpl=25;
% mu=2e9;
% nu=0.25;
% [y,z]=meshgrid(linspace(-c,c,smpl),linspace(-c,c,smpl));
% [dimx,dimy]=size(y);
% [th,r]=cart2pol(y,z);
% y(r>c)=NaN;
% z(r>c)=NaN;
% [Dn]=PennyCrackAssymGrad_disp(y(:),z(:),gamma,mu,nu,c);
% mesh(y,z,reshape(Dn,dimx,dimy),'FaceColor' ,'flat');
%
%  Author: Tim Davis     

E=(2.*mu).*(1+nu);
Dn=(4.*z.*(4./(3.*pi)).*(gamma.*c).*sqrt((c.^2 - y.^2 - z.^2)./c.^2).*(1-nu.^2))./E;
    
end