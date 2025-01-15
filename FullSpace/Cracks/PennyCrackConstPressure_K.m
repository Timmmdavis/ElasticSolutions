function [K] = PennyCrackConstPressure_K(p,c)
% PennyCrackConstPressure_K: Stress intensity around tip-line of a penny
% shaped crack subject to internal pressure.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     p             - Pressure inside the crack.
%
%     c             - Crack radius.
%
% Arguments: (output)
% 	  K             - Stress intensity around tip-line
%
%  Author: Tim Davis

%Tada P.342
K=(2/pi)*p*sqrt(pi*c);

end

