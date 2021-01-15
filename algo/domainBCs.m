function domVars=domainBCs(allVars,kLayers)
    
    domVars=endsWith(allVars,'C'); %all collar conditions return true
    isSoil=startsWith(allVars,'psiL'); %pick soil psi
    is=discoverIndices(allVars(isSoil),'psiL'); %find its layer indices
    isDom=ismember(is,kLayers); %which ones are in domain?
    domVars(isSoil)=isDom; %those return true
    
end