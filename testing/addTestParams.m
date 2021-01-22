function params=addTestParams(plant)

    params=plant.params;
    params.r=plant.R;
    params.b=plant.b;
    params.L=plant.L;
    params.kx=plant.kx;
    params.kr=plant.kr;

end