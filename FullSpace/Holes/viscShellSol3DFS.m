
function [ur,Srr,Stt]=viscShellSol3DFS(R1,R2,eta,mu,t,rIn,p,nu)
%Fullspace. 

K1=((2*mu)*(1+nu))/(3*(1-(2*nu)));
K=K1;

InShell=rIn>R1 & rIn<R2;
ur=zeros(size(rIn));
Srr=zeros(size(rIn));
Stt=zeros(size(rIn));
for i=1:numel(rIn)
    r=rIn(i);
    if InShell(i)==1  


        ur(i)=(exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(4*mu*p*r^3 + 3*K*R1^3*p - 3*K*R2^3*p - 4*R2^3*mu*p))/(12*K*mu*r^2) + (p*(3*K*R2^3 + 4*R2^3*mu - 4*mu*r^3))/(12*K*mu*r^2);
        Srr(i) =- p - (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(p*R1^3 - p*r^3))/r^3;
        Stt(i) =(exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(p*R1^3 + 2*p*r^3))/(2*r^3) - p;
 

    else

        ur(i)=(R2^3*p)/(4*mu*r^2) + (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/(4*mu*r^2);
        Srr(i)=- (R2^3*p)/r^3 - (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/r^3;
        Stt(i)=(R2^3*p)/(2*r^3) + (exp(-(3*K*R1^3*mu*t)/(3*K*R2^3*eta + 4*R2^3*eta*mu))*(p*R1^3 - p*R2^3))/(2*r^3);
    end
end
%In the fluid
ur(rIn<R1)=nan;
Srr(rIn<R1)=-p;
Stt(rIn<R1)=-p;

end