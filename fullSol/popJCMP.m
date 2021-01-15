function [C,C1,C2,C3,C4,C5]=popJCMP(parents,kr,kx,L,nL)

    a=sqrt(kr./kx);
    krL=kr.*L;
    
    b=a.*L;
    b2=(a.^2).*L;
    %th2=tanh(b/2);
    c1=sinh(b)./b;
    c2=(1-cosh(b))./b2;


    iL=(1:nL)';
    sL=nL+iL;
    
    cC=parents==0;
    nC=~cC;
    
    dKr=diag(krL);
    
    IM=diag(ones(2*nL,1));
    IM(iL,sL)=IM(iL,sL)-diag(ones(nL,1));
    IM(sub2ind(2*[nL nL],parents(nC),iL(nC)))=-1;

%% set indices of FB entries from network connectivity
    %which is better for large systems? 
    %a for loop in nL, or a matrix nL x nL?
    
%     ol=zeros(nL,1);
%     for i=1:nL
%         od=setdiff(find(parents==parents(i)),i);
%         if any(od) && nC(i)
%             ol(i)=od;
% %         else
% %             ol(i)=NaN;
%         end
%     end

    [ij,ji]=find(repmat(parents,[1 nL])==repmat(parents',[nL 1]));
    il3=ji(ij~=ji);
    ol3=ij(ij~=ji);
    ol3=ol3(nC(il3));
    il3=il3(nC(il3));
    pl3=parents(il3);
    
    il2=setdiff(find(nC),il3);
    pl2=parents(il2);
    
    
%% set the entry values
    
    %for junctions
    c3x0=(kx(pl3).*kx(il3).*c1(il3).*c2(ol3))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));
    c3x1=-(kx(il3).*(kx(pl3).*c1(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3)))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));
    c3x2=(kx(il3).*kx(ol3).*c1(il3).*c2(pl3))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));
    c3s0=(kx(pl3).*kx(il3).*c1(il3).*c2(ol3).*(c1(pl3) - 1))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));
    c3s1=-(kx(il3).*(kx(pl3).*c1(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3)).*(c1(il3) - 1))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));
    c3s2=(kx(il3).*kx(ol3).*c1(il3).*c2(pl3).*(c1(ol3) - 1))./(kx(pl3).*c1(pl3).*c2(il3).*c2(ol3) + kx(il3).*c1(il3).*c2(pl3).*c2(ol3) + kx(ol3).*c1(ol3).*c2(pl3).*c2(il3));

    %for knees
    c2x0=(kx(pl2).*kx(il2).*c1(il2))./(kx(pl2).*c1(pl2).*c2(il2) + kx(il2).*c1(il2).*c2(pl2));
    c2x1=-(kx(pl2).*kx(il2).*c1(pl2))./(kx(pl2).*c1(pl2).*c2(il2) + kx(il2).*c1(il2).*c2(pl2));
    c2s0=(kx(pl2).*kx(il2).*c1(il2).*(c1(pl2) - 1))./(kx(pl2).*c1(pl2).*c2(il2) + kx(il2).*c1(il2).*c2(pl2));
    c2s1=-(kx(pl2).*kx(il2).*c1(pl2).*(c1(il2) - 1))./(kx(pl2).*c1(pl2).*c2(il2) + kx(il2).*c1(il2).*c2(pl2));

    %for collars
    c1xc=kx(cC).*c1(cC)./c2(cC);
    c1x1=-kx(cC)./c2(cC);
    c1s1=kx(cC).*(1-c1(cC))./c2(cC);
    iC=find(cC);
    
%% Compose FB matrix
    sFB=[2*nL 2*nL+1];
    FB=zeros(sFB);
    
    %efficiency of sub2ind CAN be improved on with arithmetic operations...
    
    FB(sub2ind(sFB,iC,ones(sum(cC),1)))=c1xc;
    FB(sub2ind(sFB,iC,iC+1))=c1x1;
    FB(sub2ind(sFB,iC,nL+iC+1))=c1s1;
    
    FB(sub2ind(sFB,il2,pl2+1))=c2x0;
    FB(sub2ind(sFB,il2,il2+1))=c2x1;
    FB(sub2ind(sFB,il2,nL+pl2+1))=c2s0;
    FB(sub2ind(sFB,il2,nL+il2+1))=c2s1;
    
    FB(sub2ind(sFB,il3,pl3+1))=c3x0;
    FB(sub2ind(sFB,il3,il3+1))=c3x1;
    FB(sub2ind(sFB,il3,ol3+1))=c3x2;
    FB(sub2ind(sFB,il3,nL+pl3+1))=c3s0;
    FB(sub2ind(sFB,il3,nL+il3+1))=c3s1;
    FB(sub2ind(sFB,il3,nL+ol3+1))=c3s2;
    
    FB(sub2ind(sFB,sL,iL+1))=-krL;
    FB(sub2ind(sFB,sL,sL+1))=krL;
    
%% compose C and derivates

    C=IM*FB;

    C1=C(iL,1);
    C2=C(iL,1+iL);
    C3=C(iL,1+sL);

    C4=dKr*(eye(nL)+inv(C2)*C3);
    C5=dKr*inv(C2)*C1;


end