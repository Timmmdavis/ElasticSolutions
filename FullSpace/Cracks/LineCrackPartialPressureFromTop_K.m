function [K1Top,K1Base] = LineCrackPartialPressureFromTop_K(b,c,P,LowerTipFlag)
% LineCrackPartialPressureFromTop_K: Wall openings of a line crack line
% crack with a constant pressure. The pressure 'p' starts at the top of the
% crack and extends down to location 'b'.
% Eqs from. The stress analysis of cracks handbook. Third edition. ASME press. New York. 2000 
%
% Arguments: (input)
%     c            - Crack half-length.
%
%     P            - Constant pressure acting on the crack walls. Acts
%                    from crack's top tip (+c) down to
%                    location b.
%
%     b            - Location on the crack of the base of the pressure. Measured from the crack's centre
%
%     LowerTipFlag - If we are looking at the problem 'upside-down':
%
% Arguments: (output)
%     KITop        - Mode I stress intensity factor at the upper tip.
%
%     KIBase       - Mode I stress intensity factor at the lower tip.
%
% Diagram:
% Example 1:
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
% Example 2:
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
% Example 3: (LowerTipFlag==1)
%
%        P
%       ――――
%     K1Base
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
%    K1Top   
% 
% Example 4: (LowerTipFlag==1)
%
%        P
%       ――――
%     K1Base
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
%    K1Top   
%
%
% Example usage 1 (see diagrams above):
% c=1;
% b=0.5; 
% P=1e6;
% LowerTipFlag=0
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,P,LowerTipFlag);
%
% Example usage 2:
% c=1;
% b=-0.5; 
% P=1e6;
% LowerTipFlag=0
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,P,LowerTipFlag);
%
% Example usage 3:
% c=1;
% b=0.5; 
% P=1e6;
% LowerTipFlag=1
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,P,LowerTipFlag);
%
% Example usage 4:
% c=1;
% b=-0.5; 
% P=1e6;
% LowerTipFlag=1
% [K1Top,K1Base] = LineCrackPartialLinGradFromTop_K(b,c,P,LowerTipFlag);
%
%  Author: Tim Davis


if LowerTipFlag==1
    b=-b;
end

%Tada P.130
KIb1=1/sqrt(pi*c);
KIb2=P*c;
KIb3=pi/2-asin(b/c);
KIb4=sqrt(1-(b/c)^2);

for i=1:2

    if i==1
        K1Top=KIb1*(KIb2*(KIb3+KIb4));
    else
        K1Base=KIb1*(KIb2*(KIb3-KIb4));  
    end
    
end

if LowerTipFlag==1
    tmp=K1Top;
    K1Top=K1Base;
    K1Base=tmp;
end

end



