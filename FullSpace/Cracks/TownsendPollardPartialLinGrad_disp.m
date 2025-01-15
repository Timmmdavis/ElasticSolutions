function [Ds,Dn] = TownsendPollardPartialLinGrad_disp(b1,b2,c,gamma,z,nu,mu)
% TownsendPollardPartialLinGrad_disp: Wall openings of a line crack with
% an hourglass shaped stress gradient (linear), when across the whole crack
% this reaches gamma*c at the upper tip -gamma*c at the lower and gradient
% is 0 at the crack centre. The gradient 'gamma' starts at location b1 and
% extends down to b2. gamma is the actual gradient not the magnitude at b.
% Eq.28 of Pollard and Townsend 2017. (parts fixed)
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
%     z            - Location on the crack wall. Measured from the crack's centre
%
%     b1           - Location on the crack of the top of the linear
%                    gradient. Measured from the crack's centre
%
%     b2           - Location on the crack of the base of the linear
%                    gradient. Measured from the crack's centre
%
%     nu           - Poisson's ratio of host material.
%
%     mu           - Shear modulus of the host material.
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
% 	    gamma*c
%       ――――
% \←---|---→/ |
%  \←--|--→/  |
%      | |    |c & b1   z
%      | |b2  |			↑	
%      | |    |         | 
%      |     		     --→x 
%      |   
%      |   
%      |      
%      |   
%
% Example 2:
%
% 	    
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
%
%
% Example usage: Recreating Fig.7 Pollard and Townsend:
% c=1;
% b1=sind(45)*c;
% b2=cos(5*pi/8)*c;
% mu=1000e6;
% nu=0.25;
% gamma=1e6;
% z=linspace(-c,c,120)';
% x=zeros(size(z))
% [Ds,Dn] = TownsendPollardPartialLinGrad_disp(b1,b2,c,gamma,z,nu,mu)
% zwall_grad=[z+Ds/2;flipud(z-Ds/2)]; %merging pos & neg side
% xwall_grad=[x+Dn/2;flipud(x-Dn/2)];
% figure;
% patch(xwall_grad,zwall_grad,[0.9,0.9,0.9]);
% hold on;%grid on
% scatter(x+Dn/2,z+Ds/2,'r.');  %pos
% scatter(x-Dn/2,z-Ds/2,'b.');  %neg 
% title('Fig.7');
% ylabel('position, z (m)');
% xlabel('displacement, u_x (m)');
% xlim([-1.5e-4 1.5e-4])
% ylim([-1 1])
%
%  Author: Tim Davis


%Mapping to complex circle coords:
%Y location
VarTheta=acos(z/c); %theta, angle away from z (counter clock, complex coords) 
Z=exp(1i*VarTheta);
%Boundary conds coords. 
VarTheta1=acos(b1/c);
Z1=exp(1i*VarTheta1);
VarTheta2=acos(b2/c);
Z2=exp(1i*VarTheta2);

X=3-4*nu;
%Townsend Eq28
part1=(-gamma*c.^2)/(16*pi*1i*mu);
part2=2*(X*(Z.^-2)-(Z.^2)).*log(Z2/Z1)+(X+1)*(Z-(Z.^-1))*(Z1-(Z1.^-1)-Z2+(Z2.^-1));
%part3=(Z+(Z.^-1)).^2-(Z1-(Z1.^-1)).^2; %P&T incorrect eq
part3=(Z-(Z.^-1)).^2-(Z1-(Z1.^-1)).^2;
part4=X.*log((Z1-Z)./((Z1.^-1)-Z))-log((Z1-(Z.^-1))./((Z1.^-1)-(Z.^-1)));
%part5=(Z+(Z.^-1)).^2-(Z2-(Z2.^-1)).^2; %P&T incorrect eq
part5=(Z-(Z.^-1)).^2-(Z2-(Z2.^-1)).^2;
part6=X.*log((Z2-Z)./((Z2.^-1)-Z))-log((Z2-(Z.^-1))./((Z2.^-1)-(Z.^-1)));
uxuy=part1.*(part2+(part3.*part4)-(part5.*part6)); 

Ds=real(uxuy);
Dn=imag(uxuy);

end

