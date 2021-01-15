function i=discoverIndices(vars,str)
    hasStr=vars(startsWith(vars,str))';
    C=regexp(hasStr,'\d+','match');
    c=cat(1,C{:});
    i=str2double(c);
end