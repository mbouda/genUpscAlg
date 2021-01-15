function  vtPlant=formPlantFromVTP(vtpRoot)

    crds=reshape(vtpRoot.Points.Coordinates.data,[str2double(vtpRoot.Points.Coordinates.NumberOfComponents) vtpRoot.Dimensions.NumberOfPoints])';
    offsets=vtpRoot.Lines.offsets.data';
    ptParents=vtpRoot.Lines.connectivity.data';
    ranges=cat(2,cat(1,1,offsets(1:end-1)+1),offsets);
    nAxes=vtpRoot.Dimensions.NumberOfLines;


    for i=2:nAxes
        [~,j]=ismember(crds(ranges(i,1),:),crds,'rows');
        ptParents(ranges(i,1))=ptParents(j);
    end


    nL=vtpRoot.Dimensions.NumberOfPoints-vtpRoot.Dimensions.NumberOfLines; %+1?
    parents=zeros(nL,1);
    ax=zeros(nL,1);
    cx=zeros(nL,4);
    cy=zeros(nL,4);
    cz=zeros(nL,4);
    L=zeros(nL,1);

    count=0;
    for i=1:nAxes
        range=ranges(i,:);
            for j=range(1):(range(2)-1)
                count=count+1;
                pt0=crds(j,:);
                pt1=crds(j+1,:);
                v=pt1-pt0;
                cx(count,3)=v(1);
                cx(count,4)=pt0(1);
                cy(count,3)=v(2);
                cy(count,4)=pt0(2);
                cz(count,3)=v(3);
                cz(count,4)=pt0(3);
                L(count)=sqrt(sum(v.^2));
                ax(count)=i;
                if j==range(1)
                    ptParent=ptParents(j);
                    if ptParent>0
                        parAx=find(ranges(:,1)<=ptParent & ranges(:,2)>=ptParent);
                    else
                        parAx=1;
                    end
                    parents(count)=ptParent-(parAx-1);
                else
                    parents(count)=ptParents(j)-(i-1);
                end
            end
    end

    [M,mSys]=findMag2(parents,nL); 
    %mSys checks out: equal to nAxes

    vtPlant.nL=nL;
    vtPlant.cx=cx;
    vtPlant.cy=cy;
    vtPlant.cz=cz;
    vtPlant.parents=parents;
    vtPlant.L=L;
    vtPlant.M=M;
    vtPlant.mSys=mSys;
    vtPlant.nAxes=nAxes;
    vtPlant.AX=ax;
    vtPlant.R=vtpRoot.CellData.radius.data(ax)';

end