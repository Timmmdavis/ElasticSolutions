function [K1Top,K1Base] = TownsendPollardPartialLinGrad_K(b1,b2,c,gamma)
% TownsendPollardPartialLinGrad_disp: Stress intensity of a line crack with
% an hourglass shaped stress gradient (linear), when across the whole crack
% this reaches gamma*c at the upper tip -gamma*c at the lower and gradient
% is 0 at the crack centre. The gradient 'gamma' starts at location b1 and
% extends down to b2. gamma is the actual gradient not the magnitude at b.
% Eq.23 of Pollard and Townsend 2017.
% Fluid-filled fractures in Earth's lithosphere: Gravitational loading,
% interpenetration, and stable height of dikes and veins. Journal of
% Structural Geology. 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     gamma        - Linear stress gradient acting on the crack walls. Acts
%                    from b1 tip (+c, magnitude=gamma*b1) down to
%                    location b2 (magnitude=gamma*b2).
%
%     b1           - Location on the crack of the top of the linear
%                    gradient. Measured from the crack's centre
%
%     b2           - Location on the crack of the base of the linear
%                    gradient. Measured from the crack's centre
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%
% Diagram:
% Example 1:
%
% 	    gamma*c
%       ――――
%     Ktop
% \←---|---→/ |
%  \←--|--→/  |
%      | |    |c & b1   
%      | |b2  |				
%      | |    |          
%      |     		    
%      |   
%      |   
%      |      
%      |   
%    Kbase
%
% Example 2:
%
% 	 Ktop   
%       gamma*b1
%      | ―――   |
%  \←--|--→/ | |
%   \←-|-→/  | |c
%    \←|→/ b1| |
%     \|/    | |
%     /|\  |    
%    /→|←\ |-b2  
%   /-→|←-\|   
%      |      
%      |     
%    Kbase
%
%
% Example usage: Recreating stress intenties at the tips for Fig.7 Pollard and Townsend:
% c=1;
% b1=sind(45)*c;
% b2=cos(5*pi/8)*c;
% gamma=1e6;
% [K1Top,K1Base] = TownsendPollardPartialLinGrad_K(b1,b2,c,gamma)
%
%  Author: Tim Davis

%Mapping to complex circle coords:
%Boundary conds coords. 
VarTheta1=acos(b1/c);
Z1=exp(1i*VarTheta1);
VarTheta2=acos(b2/c);
Z2=exp(1i*VarTheta2);

%Townsend Eq2
KIa=((gamma*c)*sqrt(pi*c))/(2*pi*1i);
KIb1=log(Z2/Z1);
KIb2=(Z1-(Z1^-1)-Z2+(Z2^-1));
KIb3=0.25*((Z1^2)-(Z1^-2)-(Z2^2)+(Z2^-2)); 

K1Top=real(+KIa*(+KIb1-KIb2-KIb3));
K1Base=real(-KIa*(+KIb1+KIb2-KIb3));


end

