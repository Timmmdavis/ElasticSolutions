function [K1Top,K1Base] = TownsendPollardPartialPressure_K(b1,b2,c,P)
% TownsendPollardPartialPressure_disp: Stress intensity of a line crack
% with a constant pressure (p). The pressure starts at location b1 and
% extends down to b2.
% Eq.21 of Pollard and Townsend 2017. 
% Fluid-filled fractures in Earth's lithosphere: Gravitational loading,
% interpenetration, and stable height of dikes and veins. Journal of
% Structural Geology. 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     P            - Constant pressure acting on the crack walls. Acts
%                    from crack's top tip (+c, magnitude=gamma*c) down to
%                    location b (magnitude=gamma*b).
%
%     b1            - Location on the crack of the top of the pressure. 
%                     Measured from the crack's centre
%
%     b2            - Location on the crack of the base of the pressure. 
%                     Measured from the crack's centre
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%
% Diagram:
% Example 1:
%
% 	     P
%       ――――
%     Ktop
% |←---|---→| |
% |←---|---→| |			 
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
% 	  Ktop  
%         P
%      | ―――       |
% |←---|---→|    | |
% |←---|---→|    | |c
% |←---|---→|  b1| |
% |←---|---→|    | |
% |←---|---→| |    
% |←---|---→| |-b2  
% |←---|---→| |   
%      |      
%      |    
%    Kbase
%
%
% % Example usage: Recreating stress intenties at the tips for Pollard and Townsend Fig.6
% c=1;
% b1=sind(45)*c;
% b2=cos(5*pi/8)*c;
% P=1e6;
% [K1Top,K1Base] = TownsendPollardPartialPressure_K(b1,b2,c,P)
%
%  Author: Tim Davis

%Mapping to complex circle coords:
%Boundary conds coords. 
Theta1=acos(b1/c);
Z1=exp(1i*Theta1);
Theta2=acos(b2/c);
Z2=exp(1i*Theta2);

%Townsend Eq1
KIa=(P*sqrt(pi*c))/(2*pi*1i);
KIb1=2*log(Z2/Z1);
KIb2=(Z1-(Z1^-1)-Z2+(Z2^-1)); 

K1Top=real(KIa*(KIb1-KIb2));
K1Base=real(KIa*(KIb1+KIb2));


end

