function [K] = PennyCrackAssymGrad_K(Phi,gamma,c)
% PennyCrackAssymGrad_K: Stress intensity factor of a penny-shaped
% crack with an hourglass shaped stress gradient (linear) that reaches
% gamma*c at the upper tip -gamma*c at the lower and gradient is 0 at the
% crack centre.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack radius.
%
%     gamma        - Linear stress gradient acting on the crack walls.
%
%     Phi          - Angle away from top tip 
%
% Arguments: (output)
% 	  K            - Mode I stress intensity factor around tip (at tip-line location Phi).
%   
% Diagram 1. cross section:
%
% 	    gamma*c
%       ――――	
%    K1(Phi=0 rad)
% \←---|---→/
%  \←--|--→/
%   \←-|-→/
%    \←|→/
%     \|/
%     /|\     |
%    /→|←\    |
%   /-→|←-\   |zc (c)
%  /--→|←--\  |
% /---→|←---\ |
%     K1(Phi=pi rad)
% 	    ――――	
% 	    -gamma*c
%
% Diagram 2. Looking at fracture's face:
%
% 
%                   ooo OOO OOO ooo
%               oOO        |        OOo
%           oOO            |            OOo
%        oOO               |cz             OOo
%      oOO                 |                 OOo
%    oOO                   | Phi (Φ)          OOo
%   oOO                    |   /               OOo
%  oOO                     |__/                 OOo
% oOO                      | /                   OOo
% oOO                      |/_____________________OOo
% oOO                             cy              OOo
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
% Example usage:
% c=1;
% gamma=1000*9.81;
% Phi=linspace(0,2*pi,360);
% [K] = PennyCrackAssymGrad_K(Phi,gamma,c)
% [z,y] = pol2cart(Phi,c)
% scatter(y,z,15,K)
% colorbar
%
%  Author: Tim Davis

%Tada stress analysis of cracks handbook - P.355
K=(4/(3*pi))*gamma*c*sqrt(pi*c).*cos(Phi);

end

