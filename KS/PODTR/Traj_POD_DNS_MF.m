%%
clear;
%% Load the data.
load('MatrixFF.mat') % The FreeFem++ matrices.
load('PODMF.mat');

fid=fopen('mfcorr.txt','rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
mfcorr=b{1};
fclose(fid);

fid=fopen('data.txt','rt');
num=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
dt=fscanf(fid,'%f',1)	% time-step
step=fscanf(fid,'%d',1)	% number of time-step between two snapshots
nbre=(num-1)/step	% number of snapshots
fclose(fid);

fid=fopen(strcat('../DNS2/vect2/cbf_100001_.txt'),'rt');
num=fscanf(fid,'%d',1)
fclose(fid)

yy=zeros(size(v,2),nbre+1);

for i=1:nbre+1    % loop on snapshots of direct simulation
    i
    ii=(i-1)*step;
    fid=fopen(strcat('../DNS2/vect2/cbf_',num2str(100000+ii+1),'_.txt'),'rt');
    num=fscanf(fid,'%d',1);
    b=textscan(fid,'%f',num);
    yy(:,i)=v' * mass * (b{1}-mfcorr);
    fclose(fid);
end

%%
DNS_t=dt:(dt*step):(size(yy,2)*dt*step);
%%
save('Snap_State_POD_MF.mat','yy','DNS_t','dt','step');
