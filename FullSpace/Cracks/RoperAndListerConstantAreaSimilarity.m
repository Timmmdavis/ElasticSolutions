function [c,zpnts,hpnts,z,h]=RoperAndListerConstantAreaSimilarity(A,deltagamma,mu,nu,eta,Kc,t)
% Draw the constant area similarity solutions of Roper and Lister 2007.
% Inputs:
% A - area of crack (total including head)
% deltagamma - (rhor-rhof)g
% mu - shear mod
% nu - Young's mod
% eta - fluid visc
% Kc - fracture toughness
% t - time since start
%
% Outputs:
% c - weertman half length
% ztail - points in tail of crack
% htail - opening of one face of the tail (Dn/2)
% zhead - points in head of crack
% hhead - opening of one face of the head (Dn/2)
% 
% Example: To Draw
% figure;hold on
% plot(htail,ztail)
% plot((hhead),zhead+z+c);
% scatter(h,z,'k');
% ylabel('height from crack base')
% xlabel('width - opening D_n')
% title('Roper and Lister Fig.8 - dimensional')
% WhiteFigure


m=mu/(1-nu);

% %Checking Roper and Lister is correct on dimensional area of Weertman head
% %(not quite - wrong by constant: (pi)
% c=(Kc/(delta_gamma*sqrt(pi)))^(2/3); %Twns and Poll crit half len
% P0=0.5*delta_gamma*c;%Twns and Poll internal p
% AA=(pi*(1-nu)*P0*c^2)/(G);%Davis Healy Rivalta - Eq.3 (*2 for whole area!)
% %Area of head (weertman crack area) - Roper below Eq.6.5
% A0=((pi*Kc^2*(1-nu))/(2*G*delta_gamma))/(pi); 

A0=((Kc^2*(1-nu))/(2*mu*deltagamma)); %Area of head (weertman crack area)
at=A-A0; %Area of tail

%Half height of weertman crack: Rivalta 2015 dyke review Eq.3
c=(Kc/(deltagamma*sqrt(pi)))^(2/3);
%Eqs.6.6 to 6.8 Roper and Lister
K=((4*Kc^4*t)/(at*m^3*eta))^(1/4); %Dmlss toughness
%Eq.6.7b Roper and Lister
z=((9.*at.^2.*deltagamma.*t)./(16.*eta)).^(1/3);
%Eq.6.7a Roper and Lister
h=sqrt((eta.*z)./(deltagamma.*t)); %Opening at top of tail...

ztail=linspace(0,z,1000);
%Eq.6.7a Roper and Lister
htail=sqrt((eta.*ztail)./(deltagamma.*t));


zhead=linspace(-c,c,1000);
% %Rivalta 2015 dyke review Eq.4 % HALF OPENING!
% hhead=(((1-nu)*Kc)/(2*mu)).*sqrt(c./pi).*sqrt(1-(zhead./c).^2).*(1+(zhead./c));

p=0.5*deltagamma*c;%Twns and Poll
E=(2*mu)*(1+nu);
Eprime=E/(1-nu^2);
[DnP,~] = LineCrackConstPressure_disp(c,p,Eprime,zhead);
[DnG,~] = LineCrackAssymGrad_disp(c,-deltagamma,zhead,Eprime);
hhead=(DnP-DnG)/2;

%Put in right place relative to each other
zhead=zhead+z-c;
ztail=ztail-2*c;

%Join solutions
hhead(zhead<z-c)=[];
zhead(zhead<z-c)=[];

zpnts=[ztail,zhead]+2*c;
hpnts=[htail,hhead];

%z=z-c;

end
