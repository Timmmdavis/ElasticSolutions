
function [A]=TadaTractionIntegrationVerticalCrackA(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol)

Afun = @(b) (4*(tnFunc(b)))./(Eprime).*(sqrt(c.^2-b.^2));
A = integral(Afun,-c,c, 'RelTol', IntRelTol, 'AbsTol', IntAbsTol);

end
