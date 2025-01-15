function [KITop,KIBase] = LineCrackConstPressure_K(c,p)
% LineCrackConstPressure_K: Wall openings of a line crack 
% with a constant pressure across its faces.
% Eqs from: Tada Paris and Irwin. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     p            - Constant pressure
%
% Arguments: (output)
% 	  KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%
% Diagram:
%
% 	     p
%       ――――
%     KITop
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→|
% |←---|---→| |
% |←---|---→| |
% |←---|---→| |c
% |←---|---→| |
% |←---|---→| |
%     KIBase
%
% Example usage:
% c=1;
% p=1e6;
% [KITop,KIBase] = LineCrackConstPressure_K(c,p)
%
%  Author: Tim Davis

%Tada P.125
KITop=p*sqrt(pi*c);
KIBase=p*sqrt(pi*c);

end

