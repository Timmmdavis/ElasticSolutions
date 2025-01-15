
function [K1Top,K1Base,wf,A]=TadaTractionIntegrationVerticalCrack(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol)


[K1Top,K1Base]=TadaTractionIntegrationVerticalCrackK(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol);
[wf]=TadaTractionIntegrationVerticalCrackwf(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol);  
if nargout==4
    [A]=TadaTractionIntegrationVerticalCrackA(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol);
end

end
