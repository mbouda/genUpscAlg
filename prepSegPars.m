function [KrL,Kx,b2,c1,c2,c5]=prepSegPars(kx,kr,r,b,L)

    Kr=pi*(2*r-b).*kr./b; %no rho g: kr (opt) derived from Kr this way.
    Kx=kx.*pi.*(r-b).^2; %same here: kx derived from Lax, by dividing through with As
    a=sqrt(Kr./Kx);
    b1=a.*L;

    KrL=Kr.*L;
    b2=a.*b1;
    c1=sinh(b1)./b1;
    c2=(1-cosh(b1))./b2;
    c5=tanh(b1./2)./b1;

%     %alternatives: 
%         %Option 2 works fastest for larger systems
%           %needs nL to be added to inputs, single output: coefs
%         
% %Option 1
% %     C=cat(2,num2cell(KrL),...
% %         num2cell(b2),...
% %         num2cell(c1),...
% %         num2cell(c2),...
% %         num2cell(c5));
% %     coefs=cell2struct(C,{'KrS','b2','c1','c2','c5'},2);
%     
% %Option 2
%     coefs=struct('KrL',cell(nL,1),...
%                  'b2',cell(nL,1),...
%                  'c1',cell(nL,1),...
%                  'c2',cell(nL,1),...
%                  'c5',cell(nL,1));
%     for i=1:nL
%         coefs(i).KrS=KrL(i);
%         coefs(i).b2=b2(i);
%         coefs(i).c1=c1(i);
%         coefs(i).c2=c2(i);
%         coefs(i).c5=c5(i);
%     end
%     
% % Option 3 (incl. first code block of Opt. 2)
% %     for i=1:nL
% %         coefs(i).KrS=Kr(i).*L(i);
% %         coefs(i).b2=b2(i);
% %         coefs(i).c1=sinh(b1(i))./b1(i);
% %         coefs(i).c2=(1-cosh(b1(i)))./b2(i);
% %         coefs(i).c5=tanh(b1(i)./2)./b1(i);
% %     end
%     
    
end