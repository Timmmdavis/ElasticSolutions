function [Ds,Dn] = TownsendPollardPartialPressure_dispSPEED(b1,b2,c,P,z,nu,mu,part2,part3,part4,part5,part6)
% TownsendPollardPartialPressure_disp: Wall openings of a line crack
% with a constant pressure (P). The pressure starts at location b1 and
% extends down to b2.
% Eq.25 of Pollard and Townsend 2017. (parts fixed)
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
%     z            - Location on the crack wall. Measured from the crack's centre
%
%     b1            - Location on the crack of the top of the pressure. 
%                     Measured from the crack's centre
%
%     b2            - Location on the crack of the base of the pressure. 
%                     Measured from the crack's centre
%
%     nu            - Poisson's ratio of host material.
%
%     mu            - Shear modulus of the host material.
%
% Arguments: (output)
% 	  Dn           - The opening seperation of the fracture's walls.
%
% 	  Ds           - The shear seperation of the fracture's walls (positive
%                    is right hand wall upwards)
%
% Diagram:
% Example 1:
%
% 	     P
%       ――――
% |←---|---→| |
% |←---|---→| |			  z
%      | |    |c & b1     ↑
%      | |b2  | 	      |
%      | |    |            --→x 
%      |     
%      |   
%      |   
%      |      
%      |   
%
% Example 2:
%
% 	    
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
%
%
% % Example usage: Reproduce Pollard and Townsend Fig.6
% c=1;
% b1=sind(45)*c;
% b2=cos(5*pi/8)*c;
% mu=1000e6;
% nu=0.25;
% P=1e6;
% z=linspace(-c,c,120)';
% x=zeros(size(z))
% [Ds,Dn] = TownsendPollardPartialPressure_disp(b1,b2,c,P,z,nu,mu)
% zwall_p=[z+Ds/2;flipud(z-Ds/2)]; %merging pos & neg side
% xwall_p=[x+Dn/2;flipud(x-Dn/2)];
% figure;
% patch(xwall_p,zwall_p,[0.9,0.9,0.9]); 
% hold on;%grid on
% scatter(x+Dn/2,z+Ds/2,'r.');  %positive side
% scatter(x-Dn/2,z-Ds/2,'b.'); 
% title('Fig.6');%axis('equal')
% ylabel('position, z (m)');
% xlabel('displacement, u_x (m)');
% xlim([-8e-4 8e-4])
% ylim([-1 1])
%
%  Author: Tim Davis

% %Mapping to complex circle coords:
% %Y location
% VarTheta=acos(z/c); %theta, angle away from z (counter clock, complex coords). See Pollard and Townsend Fig.3 (here z is thier x) 
% Z=exp(1i*VarTheta);
% %Boundary conds coords. 
% Theta1=acos(b1/c);
% Z1=exp(1i*Theta1);
% Theta2=acos(b2/c);
% Z2=exp(1i*Theta2);
% 
% X=3-4*nu;
%Townsend Eq28
part1=(-P*c)/(4*pi*1i*mu);
% part2=2*(X*(Z.^-1)-Z)*log(Z2/Z1);
% part3=Z+(Z.^-1)-Z1-(Z1^-1);
% part4=X*log((Z1-Z)./((Z1^-1)-Z))-log((Z1-(Z.^-1))./((Z1^-1)-(Z.^-1)));
% part5=Z+(Z.^-1)-Z2-(Z2^-1);
% part6=X*log((Z2-Z)./((Z2^-1)-Z))-log((Z2-(Z.^-1))./((Z2^-1)-(Z.^-1)));
uxuy=part1.*(part2+(part3.*part4)-(part5.*part6));

Ds=real(uxuy);
Dn=imag(uxuy);

end

