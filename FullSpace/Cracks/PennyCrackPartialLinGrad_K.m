function [K] = PennyCrackPartialLinGrad_K(Phi,gamma,c)
% PennyCrackPartialLinGrad_K: Stress intensity factor of a penny-shaped crack with an
% hourglass shaped stress gradient (linear) that reaches gamma*c at the
% upper tip 0 at the centre. The basal part of the crack has no stress gradient.
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
%   
% Diagram 
% Example 1. cross section:
%
% 	    gamma*c
%       ――――	
%    K1(Phi=0 rad)
% \←---|---→/|
%  \←--|--→/ |
%   \←-|-→/  | zc (c)
%    \←|→/   |
%     \|/    |
%      |     
%      |     
%      |     
%      |     
%      |     
%     K1(Phi=pi rad)
% 
%
% Example 1. Looking at fracture's face (* is area where gradient acts):
%
% 
%                     OOOxOOOxooo
%               oOO********|********OOo
%           oOO************|************OOo
%        oOO***************|cz*************OOo
%      oOO*****************|*****************OOo
%    oOO*******************|*Phi (Φ)**********OOo
%   oOO********************|***/***************OOo
%  oOO*********************|__/*****************OOo
% oOO**********************|*/*******************OOo
% oOO**********************|/**********************OOo
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
%
%  Author: Tim Davis


%From z-axis 
[X,~]=pol2cart(Phi,c);
pis=pi;

K=zeros(size(Phi));
for i=1:numel(K)
    
    %Tada stress analysis of cracks handbook - summing (P.356+P.357)/2

    if X(i)>0 %Towards top of crack
        
        K(i)=gamma*c*sqrt(c*pis)*(1/(6*pis) + X(i)/c*2/(3*pis) - X(i)^2/c^2 * 4/(3*pis) +...
            4/(3*pis)*(X(i)/c + 1)^(1/2)*(X(i)/c)^(3/2));
        
    elseif X(i)<=0 %Towards base of crack
        
        K(i)=gamma*c*sqrt(c*pis)*(1/(6*pis) + X(i)/c*2/(3*pis) - X(i)^2/c^2 * 4/(3*pis) +...
            4/(3*pis)*(1 - X(i)/c)^(1/2)*(-X(i)/c)^(3/2));
        
    end
    
end

end

