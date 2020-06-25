clear
%%
load('POD.mat');
load('PODNL.mat');
load('StructFEM.mat');
load('MatrixFF.mat');
load('coord.mat');
load('StructFem.mat');

%%
q=v'*mass*v;
a=v'*lns*v;
qNL=vNL'*mass*vNL;
%% DEIM
ind = deim(vNL); % Find the DEIM indices.
sfem(ind(:))
%%
fid=fopen(strcat('deimp.txt'),'wt');
for ii=1:size(vNL,2)
    fprintf(fid,'%21.14e %21.14e %d\n',xx(ind(ii)),yy(ind(ii)),sfem(ind(ii)));
end
fclose(fid);

%%
uu=zeros(size(v));
vv=zeros(size(v));
for ii=1:size(sfem,1)
    if sfem(ii)==1
        uu(ii,:)=v(ii,:);
        uu(ii+1,:)=v(ii,:);
        vv(ii,:)=v(ii+1,:);
        vv(ii+1,:)=v(ii+1,:);
    end
end

%%

dvdx = ddx * v; % The derivatives of the POD modes.
dvdy = ddy * v; % The derivatives of the POD modes.
%%
DD=spdiags(mass2,0);
precon = @(z) z ./ DD; 
%%
qdvdx=zeros(size(v));
for i=1:size(v,2)
    qdvdx(:,i)=pcg(mass2,dvdx(:,i),1e-14,30,precon);
end
qdvdy=zeros(size(v));
for i=1:size(v,2)
    qdvdy(:,i)=pcg(mass2,dvdy(:,i),1e-14,30,precon);
end
%%
save('ROMDEIM.mat','a','ind','qdvdx','qdvdy','uu','vv');
%%
