function params=testingSetRandParams(parents)

    nL=size(parents,1);

    %here, just basic ball-park values; for realistic systems, need to set separately
    kx=repmat(5e-5,[nL 1])+1e-6*randn(nL,1);
    kr=repmat(1.5e-10,[nL 1])+1e-11*randn(nL,1); 
    b=repmat(100e-6,[nL 1])+1e-5*randn(nL,1);
    L=repmat(0.1,[nL 1])+1e-2*randn(nL,1);
    r=repmat(5e-4,[nL 1])+1e-5*randn(nL,1);

    if any(r<b) || any(b<0) || any(kx<0) || any(kr<0) || any(L<0) %check that random element did not result in unphysical setup...
        params=testingSetRandParams(parents);
    else
        [params.KrS,params.Kx,params.b2,params.c1,params.c2,params.c5]=prepSegPars(kx,kr,r,b,L);
    end
    params.nL=nL;
    params.r=r;
    params.b=b;
    params.L=L;
    params.kx=kx;
    params.kr=kr;
end