
function [wf]=TadaTractionIntegrationVerticalCrackwf(tnFunc,c,z,Eprime,IntRelTol,IntAbsTol)

%% Pressure at centre
wf=zeros(size(z));
for iii=1:numel(z)
    zi=z(iii);
    if abs(zi)>=c
        continue %No point, we know opening is zero...
    end
    wffun = @(b) (4.*(tnFunc(b)))./(pi.*Eprime).*acosh((c.^2-b.*zi)./(c.*abs(zi-b)));
    wf(iii) = integral(wffun,-c,c, 'RelTol', IntRelTol, 'AbsTol', IntAbsTol); 
end    

end
