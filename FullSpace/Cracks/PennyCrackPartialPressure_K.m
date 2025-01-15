function [K] = PennyCrackPartialPressure_K(Phi,P,c,b)
% PennyCrackPartialPressure_K: Stress intensity factor of a penny-shaped crack with an
% constant pressure that starts at the upper tip 0 extending down to location b.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack radius.
%
%     P           -  Constant pressure acting on the crack walls. Acts
%                    from crack's top tip down to
%                    location b.
%
%     Phi          - Angle away from top tip. 
%
%     b            - Base of partial pressure, measured from crack centre.
%
% Arguments: (output)
% 	  K            - Mode I stress intensity factor around tip (at tip-line location Phi).
%   
%   
% Diagram:
% Example 1 cross section:
%
%        P
%       ――――
%     K1Top
% |←---|---→| |
% |←---|---→| |
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
%    oOO                   | Phi (Φ)          OOo     |
%   oOO                    |   /               OOo    |
%  oOO                     |__/                 OOo   |b
% oOO                      | /                   OOo  |
% oOO                      |/                     OOo |
% oOO                             cy              OOo |
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
% Example 2 cross section:
%
%        P
%       ――――
%     K1Top
% |←---|---→|   |
% |←---|---→|   |
% |←---|---→|   |c
% |←---|---→|   |
% |←---|---→|   |
% |←---|---→| |    
% |←---|---→| |-b  
% |←---|---→| |   
%      |      
%      |      
%    K1Base   
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
% oOO**********************|/*********************OOo
% oOO**************************** cy *************OOo |
% oOO*********************************************OOo |
% oOO*********************************************OOo |
%  oOO*******************************************OOo  |-b  
%   oOO*****************************************OOo   |
%    oOO***************************************OOo    |
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
% Example usage 1 (see diagrams above):
% c=1;
% P=1e6;
% b=0.5
% Phi=linspace(0,2*pi,360);
% [K] = PennyCrackPartialPressure_K(Phi,P,c,b)
%
% Example usage 2 (see diagrams above):
% c=1;
% P=1e6;
% b=-0.5;
% Phi=linspace(0,2*pi,360);
% [K] = PennyCrackPartialPressure_K(Phi,P,c,b)
%
%
%  Author: Tim Davis

%Getting good sizes
if numel(b)>numel(Phi)
    Phi=ones(size(b))*Phi;
elseif numel(Phi)>numel(b)
    b=ones(size(Phi))*b;
end

if abs(b)>c
    K=zeros(size(Phi));
    return
end

%From z-axis 
[X,~]=pol2cart(Phi,c);

K=zeros(size(Phi));
for i=1:numel(K)
    
    if X(i)>b(i) %Towards top of crack
        
        K(i)=(P/(sqrt(pi*c)))*(2*c-(sqrt(c+X(i))-sqrt(X(i)-b(i)))^2);    
        
    elseif X(i)<=b(i) %Towards base of crack
        
        K(i)=(P/(sqrt(pi*c)))*((sqrt(c-X(i))-sqrt(b(i)-X(i)))^2); 
         
    end
end

% for i=1:numel(K)
%     
%     if X(i)>b(i) %Towards top of crack
%         
%         K(i)=(P/(sqrt(pi*c)))*((sqrt(c+X(i))-sqrt(X(i)-b(i)))^2);
%         
%     elseif X(i)<=b(i) %Towards base of crack
%         
%         K(i)=(P/(sqrt(pi*c)))*(2*c-(sqrt(c-X(i))-sqrt(b(i)-X(i)))^2);    
%          
%     end
% end

%% At special locations this simplifies to:
% %At top
% Ktop=(P/(sqrt(pi*c)))*((c+c)-(sqrt(2*c)-sqrt((c)-b))^2);  
% %At base
% Kbase=(P/(sqrt(pi*c)))*((c-c)+(sqrt(2*c)-sqrt(-(c*-1)+b))^2);  
% %At side:
% Kside=(P/(sqrt(pi*c)))*(2*c-(sqrt(c)-sqrt(-b))^2); 

end

