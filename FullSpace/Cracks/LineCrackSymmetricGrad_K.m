function [KITop,KIBase] = LineCrackSymmetricGrad_K(c,gamma)
% LineCrackSymmetricGrad_K: Stress intensity factor of a line crack with
% hourglass shaped stress gradient (linear) that reaches gamma*c at both
% tips and is 0 at the crack centre.
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
%    K1Top
% \←---|---→/
%  \←--|--→/
%   \←-|-→/
%    \←|→/
%     \|/
%     /|\     |
%    /←|→\    |
%   /←-|-→\   |c
%  /←--|--→\  |
% /←---|---→\ |
%     K1Base
% 	    ――――	
% 	    gamma*c
%
% Example usage:
% c=1;
% gamma=1000*9.81;
% [KITop,KIBase] = LineCrackSymmetricGrad_K(c,gamma)
%
%  Author: Tim Davis

%Due to double pos gradient Tada P.145
K1=2/pi*gamma*c*sqrt(pi*c);

KITop=K1;
KIBase=K1;

end

