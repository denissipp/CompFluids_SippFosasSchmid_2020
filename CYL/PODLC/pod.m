fid=fopen('data2.txt','rt');
num=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
dt=fscanf(fid,'%f',1)	% time-step
step=fscanf(fid,'%d',1)	% number of time-step between two snapshots
p=fscanf(fid,'%d',1)	% total number of balanced modes
nbre=(num-1)/step	% number of snapshots
fclose(fid);

fid=fopen('init.txt','rt');
in=fscanf(fid,'%d',1)	% number of iterations in direct and adjoint simulations
fclose(fid);

% determine size of snapshots
fid=fopen(strcat('../DNS2/vect2/cbf_100001_.txt'),'rt');
n=fscanf(fid,'%d',1);
fclose(fid);

% read gramian
fid=fopen('gramian.txt','rt');
C=textscan(fid,'%f');
A=reshape(C{1},nbre+1,nbre+1);
clear C;
fclose(fid);

% perform singular value decomposition
[U,S,V]=svd(A);

fid=fopen(strcat('POD/eigs.txt'),'wt');
for j=1:size(S,1)	% loop on snapshots of direct simulation
    fprintf(fid,'%21.14e\n',S(j,j));
end
fclose(fid)

for j=1:p	% save direct balanced modes
    j
    fid=fopen(strcat('PODTIME/mode_',num2str(j),'_.txt'),'wt');
    fprintf(fid,'%d\n',nbre+1);
    fprintf(fid,'%21.14e\n',sqrt(S(j,j))*V(:,j));
    fclose(fid);
end

% compute balanced modes
v=zeros(n,p);
for i=0:nbre	% loop on snapshots of direct simulation
    i
    ii=i*step;
    fid=fopen(strcat('../DNS2/vect2/cbf_',num2str(100000+in+ii+1),'_.txt'),'rt');
    num=fscanf(fid,'%d',1);
    b=textscan(fid,'%f',num);
    fclose(fid);
    v=v+b{1}*V(i+1,1:p)*S(1:p,1:p)^(-1/2);	% update balanced mode
end
    
for j=1:p	% save direct balanced modes
    j
    fid=fopen(strcat('POD/mode_',num2str(j),'_.txt'),'wt');
    fprintf(fid,'%d\n',n);
    fprintf(fid,'%21.14e\n',v(:,j));
    fclose(fid);
end

save('POD.mat', 'v', '-v7.3');
