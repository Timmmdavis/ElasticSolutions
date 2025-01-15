function [KITop,KIBase] = LineCrackSymmetricGradQuadratic_K(c,gamma)
% LineCrackSymmetricGradQuadratic_K: Stress intensity factor of a line crack with
% hourglass shaped stress gradient (quadratric) that reaches gamma at both
% tips and is 0 at the crack centre.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     gamma        - Quadratic stress gradient value at tips. Increases
%                    along crack by: sxx=gamma*(abs(z)/a)^2; Tada p.148.
%     
%     
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%     
% Diagram:
%
%    K1Top
%      ----
% \←--|--→/
%    ←|→
%   ⎞←|→⎛
%    |||     |
%    ←|→     |
%   ⎠←|→⎝   |c
% /←--|--→\  |
%
%     K1Base
% 	    ――――	
% 	    gamma
%
% Example usage:
% c=1;
% gamma=1000*9.81;
% [KITop,KIBase] = LineCrackSymmetricGrad_K(c,gamma)
%
%  Author: Tim Davis

%Due to double pos gradient Tada P.148 
K1=(gamma*sqrt(c*pi))/2;

KITop=K1;
KIBase=K1;

end

