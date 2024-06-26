%% startup DGDD-PFVS
addpath(genpath(pwd))
%%

load geopack
domain=geopack;
%%

domain=readTiffStackYDW('C:\Users\Whydee\Downloads\saturation_pc_278.7133858.tif');
domain=imresize3(domain,0.5,'nearest');
domain(domain>0)=1;
voxelSize=1e-6;
alpha=1;
gradk=0;
FDGPA=1;
FDGPAYD=0;
micropore=0;
Pin=1;
Pout=0;
condFlag=1;
cond=double(~domain);
[K,vel,P,phi,phiEff]= solvePFVS(domain,voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond);