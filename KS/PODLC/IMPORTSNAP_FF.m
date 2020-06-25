clear all;close all;clc;

%% Read measurements of DNS Step=1
fid=fopen('../DNS2/mes_0.txt','rt');
CC=textscan(fid,'%f');
M=reshape(CC{1},62,size(CC{1},1)/62);
DNS_t2=M(1,:);
mesDNSbrut=M(2:62,:);
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
fid=fopen(strcat('x.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
xx=b{1};
fclose(fid);
save('coord.mat','xx');
%%
%%
fid=fopen(strcat('offset.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
offset=b{1};
fclose(fid);
save('offset.mat','offset');
