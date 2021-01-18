function c=discoverIndices(vars,str)
    hasStr=vars(startsWith(vars,str))';
    
    nHas=length(hasStr);
    c=zeros(nHas,1);
    for i=1:nHas
        c(i)=str2double(hasStr{i}(length(str)+1:end));
    end
    
end