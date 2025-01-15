function [K1Top,K1Base] = LineCrackPartialLinGradFromTop_KZeroAtb(b,c,gamma,LowerTipFlag)
% LineCrackPartialLinGradFromTop_K: Wall openings of a line crack line crack
% with an hourglass shaped stress gradient (linear), when across the whole
% crack this reaches gamma*c at the upper tip -gamma*c at the lower and
% gradient is 0 at the crack centre. The gradient 'gamma' starts at the top
% of the crack and extends down to 'b', gamma is the gradient not the
% magnitude of force at b.
% Eqs simplfied by Tim Davis. Simplifed form of Eq.21 of Pollard and Townsend 2017.
% Fluid-filled fractures in Earth's lithosphere: Gravitational loading,
% interpenetration, and stable height of dikes and veins. Journal of
% Structural Geology. 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     gamma        - Linear stress gradient acting on the crack walls. Acts
%                    from crack's top tip (+c, magnitude=gamma*b) down to
%                    location b (magnitude=0).
%
%     b            - Location on the crack of the base of the linear
%                    gradient. Measured from the crack's centre
%
%     LowerTipFlag - If we are looking at the problem 'upside-down':
%                    Note when this is 1, the sign of gamma*z stays the
%                    same (negative below crack centre)
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%
% Diagram:
% Example 1:
%
%       gamma*(c-b)
%       ――――
%     K1Top
%    \←|→/    |
%     \|/     |
%      | |    |c
%      | |b   |
%      | |    |
%      |     
%      |   
%      |   
%      |      
%      |  
%    K1Base   
%
% Example 2:
%
% 	    gamma*(c-b)
%       ――――
%     K1Top
% \←---|---→/|
%  \←--|--→/ |
%   \←-|-→/  |c
%    \←|→/   |
%     \|/  |-b  
%      |      
%      |   
%      |        
%    K1Base   
%
% Example 3:
%
%       gamma*b
%       ――――
%     K1Base
%    \→|←/    |
%     \|/     |
%      | |    |c
%      | |b   |
%      | |    |
%      |     
%      |   
%      |   
%      |      
%      |  
%    K1Top   
%
% Example 4:
%
% 	    gamma*c
%       ――――
%     K1Base
% \→---|---←/|
%  \→--|--←/ |
%   \→-|-←/  |c
%    \→|←/   |
%     \|/  |-b  
%      |      
%      |   
%      |        
%    K1Top   
%
% Example usage 1 (see diagrams above):
% c=1;
% b=0.5; 
% gamma=1000*9.81;
% LowerTipFlag=0
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,gamma,LowerTipFlag);
%
% Example usage 2:
% c=1;
% b=-0.5; 
% gamma=1000*9.81;
% LowerTipFlag=0
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,gamma,LowerTipFlag);
%
% Example usage 3:
% c=1;
% b=0.5; 
% gamma=1000*9.81;
% LowerTipFlag=1
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,gamma,LowerTipFlag);
%
% Example usage 4:
% c=1;
% b=-0.5; 
% gamma=1000*9.81;
% LowerTipFlag=1
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,gamma,LowerTipFlag);
%
%  Author: Tim Davis

if LowerTipFlag==1
    b=-b;
end

for i=1:2
    
    if i==1
        K1Top=(c^(1/2)*gamma)/(2*pi^(1/2))*((2*c-b)*(( c^2- b^2 )/c^2)^(1/2) - b*pi + 2*b*asin(b/c) + c*acos(b/c));
    else
        K1Base=(c^(1/2)*gamma)/(2*pi^(1/2))*((2*c+b)*((c^2- b^2 )/c^2)^(1/2) - b*pi + 2*b*asin(b/c) - c*acos(b/c));
    end 

end

%Note we assume that its positive gamma going from the base 'down' to the
%upper tip...
if LowerTipFlag==1
    tmp=-K1Top;
    K1Top=-K1Base;
    K1Base=tmp;
end

end



