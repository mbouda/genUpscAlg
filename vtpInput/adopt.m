function newPlant=adopt(plant,dtr,par)

    newPlant=plant;
    
    x0=sum(plant.cx(par,:),2);
    x1=sum(plant.cx(dtr,:),2);
    newPlant.cx(dtr,:)=[x1-x0 x0];
    
    y0=sum(plant.cy(par,:),2);
    y1=sum(plant.cy(dtr,:),2);
    newPlant.cy(dtr,:)=[y1-y0 y0];
    
    z0=sum(plant.cz(par,:),2);
    z1=sum(plant.cz(dtr,:),2);
    newPlant.cz(dtr,:)=[z1-z0 z0];
    
    newPlant.L(dtr,:)=sqrt(sum([newPlant.cx(dtr,1) newPlant.cy(dtr,1) newPlant.cz(dtr,1)].^2));
    
    newPlant.parents(dtr)=par;

%     %actually, magnitude can just be re-read from parents...
%        %then transfer daughters to it
%        %and that code will also have to adjust the M
%     if plant.parents(dtr)==0
%        plant.mSys
%        plant.M(dtr)
%         %magnitude of the offshoots does not change
%         %magnitude along the main stem must decrease by that number at each junction
%         
%     else
%         
%         
%     end
    
    
end