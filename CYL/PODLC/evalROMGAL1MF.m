clear;
%%
load('Snap_State_POD_MF.mat');
%%
load('PODMF.mat');
load('ROMDEIMMF.mat');
load('MatrixFF.mat');
%%
nROM=14
%%
v_ROM=v(:,1:nROM);
osROM=os(1:nROM);
aROM=a(1:nROM,1:nROM);
ev=eig(aROM);
[~,ii]=max(real(ev));
uev=real(ev(ii))+1i*abs(imag(ev(ii)));
tev=0.006531009709731875+1i*1.05294062310488
erroreig=abs(uev-tev)/abs(tev)*100
%%
fROMi=zeros(nROM,nROM,nROM);
for j=1:nROM
    for k=1:nROM
        vect=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j)+vv(:,j).*qdvdy(:,k)+vv(:,k).*qdvdy(:,j));
    	fROMi(j,k,:)=v_ROM'*(mass*vect);   
    end
end
%% Cure
%for i=1:nROM
%    fROMi(i,i,i)=0;
%end
%%
fROMa=zeros(nROM,nROM,nROM);
for j=1:nROM
    for k=1:nROM
		for i=1:nROM
			fROMa(j,k,i)=(5*fROMi(j,k,i)-fROMi(j,i,k)-fROMi(i,j,k)-fROMi(i,k,j)-fROMi(k,i,j)-fROMi(k,j,i))/6;
        end
    end
end
%%
fROMs=fROMi-fROMa;
num=0; den=0;
for i=1:nROM
    num=num+norm(fROMs(:,:,i),'fro')^2;
    den=den+norm(fROMi(:,:,i),'fro')^2;
end
assymetry=sqrt(num/den)*100
%%
global fROM;
in=15000
%in=7500
THor=15001
fROM=fROMi; disp 'fROMi'
%fROM=fROMa; disp 'fROMa'

rangeDNSbrut=(in+1):(in+THor);
t0=rangeDNSbrut*dt;
model_t = t0;

z = zeros(nROM, length(model_t)); % The solution array.
z(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition.

%Use 1st-order accurate time scheme for the first step.
%for j = 2
%z(:, 2) = (eye(nROM) - dt * aROM) \ (z(:, 1) + dt * dydt_nonlinear(z(:, 1)) + dt * osROM);
%end
for j = 2
z(:, 2) = (eye(nROM) - dt * aROM) \ (z(:, 1));
end

%Use 2nd-order accurate time scheme for all subsequent steps.
%for j = 3:length(model_t)
%    z(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
%        (4 * z(:, j-1) - z(:, j-2) ...
%        + 4 * dt * dydt_nonlinear(z(:, j-1))- 2 * dt * dydt_nonlinear(z(:, j-2)) + 2 * dt * osROM);
%end
for j = 3:length(model_t)
    z(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * z(:, j-1) - z(:, j-2));
end

%%
figure(2)
plot(DNS_t,yy(1,:),'r')
hold on
plot(model_t,z(1,:),'k')
hold off
%%
% Error on full snapshot
rangeDNS=in/step+1:(in/step+THor/step+1);
rangeMODEL=1:step:THor;
yye=yy(:,rangeDNS);
ye=z(:,rangeMODEL);
A=norm(ye(1:nROM,:)-yye(1:nROM,:),'fro')^2/size(yye,2); % mean error per snapshot
B=norm(yye(nROM+1:end,:),'fro')^2/size(yye,2);
CC=norm(yye(:,:),'fro')^2/size(yye,2);
errMODEL=sqrt(A/CC)*100
errTRUNC=sqrt(B/CC)*100
