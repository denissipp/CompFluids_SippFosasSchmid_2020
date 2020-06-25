clear
%%
load('PODMF.mat');
load('PODNLMF.mat');
load('MatrixFF.mat');
load('coord.mat');
load('StructFem.mat');

%%
q=v'*mass*v;
diag(q)
a=v'*lnsmf*v;
eig(a)
qNL=vNL'*mass*vNL;
%% DEIM
ind = deim(vNL); % Find the DEIM indices.
sfem(ind(:))
%%
fid=fopen(strcat('deimpMF.txt'),'wt');
for ii=1:size(vNL,2)
    fprintf(fid,'%21.14e %d\n',xx(ind(ii)),sfem(ind(ii)));
end
fclose(fid);

%%
uu=zeros(size(v));
for ii=1:size(sfem,1)
    if sfem(ii)==1
        uu(ii,:)=v(ii,:);
    end
end

%%

dvdx = ddx * v; % The derivatives of the POD modes.
%%
DD=spdiags(mass2,0);
precon = @(z) z ./ DD; 
%%
qdvdx=zeros(size(v));
for i=1:size(v,2)
    qdvdx(:,i)=pcg(mass2,dvdx(:,i),1e-14,30,precon);
end
%%
load('offset.mat');
os=v'*offset;
%%
save('ROMDEIMMF.mat','a','os','ind','qdvdx','uu');
%%
