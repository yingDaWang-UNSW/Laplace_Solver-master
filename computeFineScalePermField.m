%% Load RAW file
fname = "C:\Users\yingd\Downloads\geopack 3.raw";
dtype = 'uint8';
[Nx, Ny, Nz] = deal(650, 650, 1000);
fileID = fopen(fname, 'r');
data = fread(fileID, Nx*Ny*Nz, [dtype, '=>', dtype]);
fclose(fileID);
data = reshape(data, [Nx, Ny, Nz]);

load geopack;
data=geopack;
[Nx, Ny, Nz] = deal(150, 150, 150);

%% Define block size
block_size = 50;

% Define how to split each dimension
x_sizes = repmat(block_size, 1, floor(Nx / block_size));
if mod(Nx, block_size) ~= 0
    x_sizes(end + 1) = mod(Nx, block_size);
end

y_sizes = repmat(block_size, 1, floor(Ny / block_size));
if mod(Ny, block_size) ~= 0
    y_sizes(end + 1) = mod(Ny, block_size);
end

z_sizes = repmat(block_size, 1, floor(Nz / block_size));
if mod(Nz, block_size) ~= 0
    z_sizes(end + 1) = mod(Nz, block_size);
end

%% Split into sub-blocks
data_blocks = mat2cell(data, x_sizes, y_sizes, z_sizes);

%% run pfvs for each subblock 3 times
voxelSize=1e-5;
alpha=1; % for porous media, this should be 1, for fractured media, it should be 2
gradk=0;
FDGPA=1;
FDGPAYD=0;
micropore=0;
Pin=1;
Pout=0;
condFlag=0;
cond=[];
lx=size(data_blocks,1);
ly=size(data_blocks,2);
lz=size(data_blocks,3);

permField=cell(lx, ly, lz);

for i=1:lx
    for j=1:ly
        for k=1:lz
            [Kz,vel,P,phi,phiEff]= solvePFVS(data_blocks{i,j,k},voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond);
            [Kx,vel,P,phi,phiEff]= solvePFVS(permute(data_blocks{i,j,k},[2,3,1]),voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond);
            [Ky,vel,P,phi,phiEff]= solvePFVS(permute(data_blocks{i,j,k},[1,3,2]),voxelSize,alpha,gradk,FDGPA,FDGPAYD,micropore,Pin,Pout,condFlag,cond);
            permField{i,j,k}=[Kx(1), Ky(1), Kz(1)];
        end
    end
end


%% run again to get overall perm

cond=[];

for i=1:lx
    for j=1:ly
        for k=1:lz
            cond=[cond;permField{i,j,k}];
        end
    end
end
[K,vel,P,phi,phiEff]= solvePFVS(ones(lx, ly, lz),1,0,0,0,0,1,1,0,1,cond.*1e-12);


