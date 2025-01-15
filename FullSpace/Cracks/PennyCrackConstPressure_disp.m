function [Dn]=PennyCrackConstPressure_disp(y,z,p,mu,nu,c)
% PennyCrackAssymGradient_disp: Compute the separation of the cracks faces
% (D_n) due to internal pressure p.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000
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
%     p             - Internal pressure acting on faces.
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
% 2. Constant pressure boundary condition, 2D cross section (xz):
%
% 	   
% 	     p
%       ――――
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→| |
% |←---|---→| |
% |←---|---→| |c  (z-axis)
% |←---|---→| |
% |←---|---→| |
%
% p=1e6;
% c=500;
% smpl=25;
% mu=2e9;
% nu=0.25;
% [y,z]=meshgrid(linspace(-c,c,smpl),linspace(-c,c,smpl));
% [dimx,dimy]=size(y);
% [th,r]=cart2pol(y,z);
% y(r>c)=NaN;
% z(r>c)=NaN;
% [Dn]=PennyCrackConstPressure_disp(y,z,p,mu,nu,c)
% mesh(y,z,reshape(Dn,dimx,dimy),'FaceColor' ,'flat');
%
%  Author: Tim Davis    

E=(2.*mu).*(1+nu);

%Tada 342 - Opening due to pressure
r=(y.^2+z.^2).^(0.5);
Dn=((8.*(1-nu.^2))./(pi.*E)).*p.*sqrt(c.^2-r.^2);

end