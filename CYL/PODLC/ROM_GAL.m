clear;
%%
load('MatrixFF.mat') % The FreeFem++ matrices.
%%
load('POD.mat');
%%
q=v'*mass*v;
a=v'*lns*v;
%% Construct the components of the reduced-order model.
eig(a)
figure(1)
plot(eig(a),'x');
%%
%%
dvdx = ddx * v; % The derivatives of the POD modes.
dvdy = ddy * v; % The derivatives of the POD modes.
%%
load('StructFem.mat');
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
f=zeros(size(v,2),size(v,2),size(v,2));
for i=1:size(v,2)
    i
    for j=1:size(v,2)
        for k=1:size(v,2)
            ww=v(:,i).*(uu(:,j).*dvdx(:,k)+vv(:,j).*dvdy(:,k));
            f(j,k,i)=sum(ww,1);
        end
    end
end
%%
save('ROMGAL.mat','a','f');
%%
