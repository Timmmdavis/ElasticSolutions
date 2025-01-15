function [Srr,Stt,Srt,ur] = KirschSolution(p,a,r,mu)
%Kirsch pressure only - P.248 my old Jaeger and Cook
Srr=p.*a.^2./r.^2;
Stt =-p.*a.^2./r.^2;
Srt=0;
ur=(p.*a.^2)./(2.*mu.*r);

end