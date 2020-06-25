clear;
%%
load('MatrixFF.mat') % The FreeFem++ matrices.
%%
load('PODMF.mat');
%%
q=v'*mass*v;
a=v'*lnsmf*v;
%% Construct the components of the reduced-order model.
eig(a)
figure(1)
plot(eig(a),'x');
%%
%%
dvdx = ddx * v; % The derivatives of the POD modes.
%%
load('StructFem.mat');
%%
uu=zeros(size(v));
for ii=1:size(sfem,1)
    if sfem(ii)==1
        uu(ii,:)=v(ii,:);
    end
end

%%
f=zeros(size(v,2),size(v,2),size(v,2));
for i=1:size(v,2)
    i
    for j=1:size(v,2)
        for k=1:size(v,2)
            ww=v(:,i).*(uu(:,j).*dvdx(:,k));
            f(j,k,i)=sum(ww,1);
        end
    end
end
%%
load('offset.mat');
os=v'*offset;
%%
save('ROMGALMF.mat','a','f','os');
%%
