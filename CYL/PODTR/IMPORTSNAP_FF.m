clear all;close all;clc;
%%
%% Read measurements of DNS Step=1
fid=fopen('../DNS2/mes_0.txt','rt');
CC=textscan(fid,'%f');
M=reshape(CC{1},101,size(CC{1},1)/101);
DNS_t2=M(1,:);
mesDNSbrut=M(2:101,:);
clear CC,M;
save('mes_0.mat','DNS_t2','mesDNSbrut');
fclose(fid);%%
%%
fid=fopen(strcat('structfem.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
sfem=b{1};
fclose(fid);
save('StructFem.mat','sfem');
%%
fid=fopen(strcat('offset.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
offset=b{1};
fclose(fid);
save('offset.mat','offset');
%%
%%
fid=fopen(strcat('mfcorr.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
mfcorr=b{1};
fclose(fid);
save('mfcorr.mat','mfcorr');
%%
%%
fid=fopen(strcat('x.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
xx=b{1};
fclose(fid);
fid=fopen(strcat('y.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
yy=b{1};
fclose(fid);
save('coord.mat','xx','yy');
%%
unstable=1;
fid=fopen(strcat('../Eigs/ev0.txt'),'rt');
num=fscanf(fid,'%d',1)
fclose(fid)

vu=zeros(num,2*unstable);

for i=1:unstable    % loop on snapshots of direct simulation
    i
    fid2=fopen(strcat('../Eigs/ev',num2str(i-1),'.txt'),'rt');
    num=fscanf(fid2,'%d',1);
    b=textscan(fid2,'(%f,%f)',num);
    vu(:,2*(i-1)+1)=b{1};
    vu(:,2*(i-1)+2)=b{2};
    fclose(fid2);
end

wu=zeros(num,2*unstable);

for i=1:unstable    % loop on snapshots of direct simulation
    i
    fid2=fopen(strcat('../Eigs/ea_scaled',num2str(i-1),'.txt'),'rt');
    num=fscanf(fid2,'%d',1);
    b=textscan(fid2,'(%f,%f)',num);
    wu(:,2*(i-1)+1)=2.*b{1};
    wu(:,2*(i-1)+2)=2.*b{2};
    fclose(fid2);
end
save('UnstableBasis.mat','vu','wu');
