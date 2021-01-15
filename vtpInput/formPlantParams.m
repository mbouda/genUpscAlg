function params=formPlantParams(plant)

    [params.KrS,params.Kx,params.b2,params.c1,params.c2,params.c5]=...
        prepSegPars(plant.kx,plant.kr,plant.R,plant.b,plant.L);
    params.nL=plant.nL; %currently duplicate, but this is where the algorithm will look for it
    
    
end