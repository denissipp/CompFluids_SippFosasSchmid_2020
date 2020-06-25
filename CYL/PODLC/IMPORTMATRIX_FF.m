clear all;close all;clc;

%%
%IMPORT MATRICES
%Import LNS matrix
fid=fopen('LNS.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
lns=sparse(C{1},C{2},C{3});
clear C;
%%
%Import LNSMF matrix
fid=fopen('LNSMF.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
lnsmf=sparse(C{1},C{2},C{3});
clear C;
%%
%IMPORT MASS Matrix
fid=fopen('M.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
mass=sparse(C{1},C{2},C{3});
clear C;
%%
%IMPORT MASS Matrix
fid=fopen('M2.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
mass2=sparse(C{1},C{2},C{3});
clear C;
%%
%IMPORT DErivative matrix
fid=fopen('DerX.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
ddx=sparse(C{1},C{2},C{3});
clear C;
%%
%IMPORT DErivative matrix
fid=fopen('DerY.txt','rt');
line = fgetl(fid);
line = fgetl(fid);
line = fgetl(fid);
num=fscanf(fid,'%d %d %d %d',4)
C=textscan(fid,'%f %f %f');
ddy=sparse(C{1},C{2},C{3});
clear C;
%%
clear line fid num

save('MatrixFF.mat');
