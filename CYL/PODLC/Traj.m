clear;
%% Load the data.
load('MatrixFF.mat') % The FreeFem++ matrices.
load('POD.mat');
%%
fid=fopen(strcat('mfcorr.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
mfcorr=-b{1};
fclose(fid);
%%
fid=fopen(strcat('../DNS/vect/cbf_100001_.txt'),'rt');
num=fscanf(fid,'%d',1)
fclose(fid)
%%
eigen=zeros(num,1);

fid=fopen(strcat('../Eigs/evr1.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
eigen=b{1};
fclose(fid);
%%
fid=fopen(strcat('../Eigs/evi1.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
eigen = eigen+b{1}*1j;
fclose(fid);
%%
eigenMF=zeros(num,1);

fid=fopen(strcat('../EigsMF/evr0.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
eigenMF=b{1};
fclose(fid);
%%
fid=fopen(strcat('../EigsMF/evi0.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
eigenMF = eigenMF+b{1}*1j;
fclose(fid);
%%
save('SnapBF.mat','mfcorr','eigen','eigenMF');
