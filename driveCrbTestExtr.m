addpath ./algo ./fullSol ./testing ./vtpInput

kx=5e-5;
kr=1.5e-10;
b=100e-6;
nLayInit=8;

dataDir='/run/media/mbouda/OS/Users/Martin/Documents/graphics/';

rsaFile='Lupinus_angustifolius_Chen_2011_LAB201008_7denni_simulace.vtp';
[plant,zMin,zLims,dz]=importPlant(strcat(dataDir,rsaFile),nLayInit,kr,kx,b); 
plant.params=addTestParams(plant);


testSet(1).parents=plant.parents;
testSet(1).inLayer=plant.inLayer;
testSet(1).params=plant.params;

dataDir='./testing/crbTestSet/';

nDay=[10:3:16 22:3:28 34:3:40];
for i=1:size(nDay,2)
    rsaFile=sprintf('RLab_210117_Lupinus_angustifolius_Chen_2011_%ddenni_simulace.vtp',nDay(i));
    [plant,zMin,zLims,dz]=importPlant(strcat(dataDir,rsaFile),nLayInit,kr,kx,b); 
    plant.params=addTestParams(plant);
    testSet(i+1).parents=plant.parents;
    testSet(i+1).inLayer=plant.inLayer;
    testSet(i+1).params=plant.params;
end

dataDir='/run/media/mbouda/OS/Users/Martin/Documents/graphics/';
rsaFile='Lupinus_angustifolius_Chen_2011_LAB201008_42denni_simulace.vtp';
[plant,zMin,zLims,dz]=importPlant(strcat(dataDir,rsaFile),nLayInit,kr,kx,b); 
plant.params=addTestParams(plant);
testSet(end+1).parents=plant.parents;
testSet(end).inLayer=plant.inLayer;
testSet(end).params=plant.params;
    
save('./testing/crbTestSet/testSet','testSet');

%rsaFile=sprintf('RLab_210119_Pisum_sativum_a_Pagès_2014_%ddenni_simulace.vtp',nDay);
%rsaFile='RLab_210115_workshop.vtp';
