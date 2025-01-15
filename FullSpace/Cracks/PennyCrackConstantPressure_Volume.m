function [Volume] = PennyCrackConstantPressure_Volume(nu,E,pressure,c)
% PennyCrackConstantPressure_Volume: Volume of a penny shaped crack subject
% to internal pressure
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     nu            - Poisson's ratio of host material.
%
%     E             - Young's modulus of the host material.
%
%     pressure      - Pressure inside the crack.
%
%     c             - Crack radius.
%
% Arguments: (output)
% 	  Volume        - Volume of the opening between the crack faces ('void space').
%
%  Author: Tim Davis

%Tada P.342
Volume=((16*(1-nu^2))/(3*E))*pressure*c^3;

end

