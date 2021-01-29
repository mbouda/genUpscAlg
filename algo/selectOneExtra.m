function keepExtra=selectOneExtra(eqs)

    allCs=cat(1,{eqs(:).coefs});
    nCs=cellfun(@(x)size(x,2),allCs);
    hasMax=nCs==max(nCs);
    
    minCs=cellfun(@(x)min(abs(x)),allCs)';
    [~,iMax]=min(minCs(hasMax));  %switched function from max to min
    
    iEq=find(hasMax);
    
    keepExtra=false(size(eqs,1),1);
    keepExtra(iEq(iMax))=true;
    
end