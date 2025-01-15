function [Volume] = PennyCrackPartialPressure_Volume(nu,E,P,c,b)
% PennyCrackPartialPressure_Volume: Volume of a penny shaped crack subject
%to a partial pressure from upper tip down to location b (b measured from
%centre of crack).
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     nu            - Poisson's ratio of host material.
%
%     E             - Young's modulus of the host material.
%
%     c            - Crack radius.
%
%     P           -  Constant pressure acting on the crack walls. Acts
%                    from crack's top tip (+c, magnitude=gamma*c) down to
%                    location b (magnitude=gamma*b).
%
%     b            - Base of partial pressure, measured from crack centre.
%
% Arguments: (output)
% 	  Volume        - Volume of the opening between the crack faces ('void space').
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
% nu=0.25;
% E=1e9;
% [Volume] = PennyCrackPartialPressure_Volume(nu,E,P,c,b)
%
% Example usage 2 (see diagrams above):
% c=1;
% P=1e6;
% b=-0.5;
% nu=0.25;
% E=1e9;
% [Volume] = PennyCrackPartialPressure_Volume(nu,E,P,c,b)
%
%
%  Author: Tim Davis

%Tada P.354
Volume=((4*(1-nu^2))/(3*E))*P*(c-b)^2*(2*c+b);

end

