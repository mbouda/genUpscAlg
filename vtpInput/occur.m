function [o,nodeType]=occur(A)
%o=occur(A)
%count occurances of non-negative integer entries in column vector A

    A=A+1; %to get rid of 0, with no effect on outcome
    szA=size(A,1);
    mxA=max(A);
    
    
    B=sparse(A,1:szA,true,mxA,szA);
    C=sum(B,2);
    D=find(C>0);
    nodeType=[D-1 full(C(D))];
    nodeType(1,:)=[];
    o=full(C(A));
end