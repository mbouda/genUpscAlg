function plant=adjustVTPlant2(plant)


    plant.cx=plant.cx/100; %these are clearly in cm coming from CPlantBox
    plant.cy=plant.cy/100;
    plant.cz=plant.cz/100;
    %they are also actually linear!
    
    plant.cx=plant.cx(:,3:4);
    plant.cy=plant.cy(:,3:4);
    plant.cz=plant.cz(:,3:4);
    
    plant.L=plant.L/100;
    plant.R=plant.R/100;
end