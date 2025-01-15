function [KITop,KIBase] = LineCrackAssymGrad_K(c,gamma)
% LineCrackAssymGrad_K: Stress intensity factor of a line crack with an
% hourglass shaped stress gradient (linear) that reaches gamma*c at the
% upper tip -gamma*c at the lower and gradient is 0 at the crack centre.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     gamma        - Linear stress gradient acting on the crack walls.
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%   
% Diagram:
%
% 	   
% 	    gamma*c
%       ――――	
%    K1Top
% \←---|---→/
%  \←--|--→/
%   \←-|-→/
%    \←|→/
%     \|/
%     /|\     |
%    /→|←\    |
%   /-→|←-\   |c
%  /--→|←--\  |
% /---→|←---\ |
%     K1Base
% 	    ――――	
% 	    -gamma*c
%
% Example usage:
% c=1;
% gamma=1000*9.81;
% [KITop,KIBase] = LineCrackAssymGrad_K(c,gamma)
%
%  Author: Tim Davis

%Tada P.153
KITop=0.5*(gamma*c)*sqrt(pi*c);
KIBase=-0.5*(gamma*c)*sqrt(pi*c);

end

