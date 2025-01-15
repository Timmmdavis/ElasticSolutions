
function [K1Top,K1Base]=TadaTractionIntegrationVerticalCrackK(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol)

Kfun = @(b) (tnFunc(b))./(sqrt(pi.*c)).*(sqrt(c.^2-b.^2)./(c-b));
K1Top = integral(Kfun,-c,c, 'RelTol', IntRelTol, 'AbsTol', IntAbsTol);

% Kfun = @(b) (TnFunTotal(b))./(sqrt(pi.*c)).*(sqrt(c.^2-b.^2)./(c+b));
% K1Base = integral(Kfun,-c,c, 'RelTol', IntRelTol, 'AbsTol', IntAbsTol);
K1Base=K1Top; %Assuming symmetry


end
