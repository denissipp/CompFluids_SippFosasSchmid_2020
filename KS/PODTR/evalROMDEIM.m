clear;
%%
load('Snap_State_POD.mat');
%%
load('POD.mat');
load('PODNL.mat');
load('ROMDEIM.mat');
load('MatrixFF.mat');
%%
load('mes_0.mat')
%% Read measurements of POD
fid=fopen('mespod.txt','rt');
CC=textscan(fid,'%f');
C=reshape(CC{1},61,size(a,1));
clear CC;
fclose(fid);%%
%%
nROM=60
nROM2=70
%%
v_ROM=v(:,1:nROM);
aROM=a(1:nROM,1:nROM);
ev=eig(aROM);
[~,ii]=max(real(ev));
uev=real(ev(ii))+1i*abs(imag(ev(ii)));
tev=0.3378487570072862+1i*0.618195700937916;
erroreig=abs(uev-tev)/abs(tev)*100
%%
vNL_ROM=vNL(:,1:nROM2);
indROM = ind(1:nROM2);
MatInv=pinv(vNL_ROM(indROM, :));
%%
fROMi=zeros(nROM,nROM,nROM);
for j=1:nROM
    for k=1:nROM
		nl=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j));
        nlROM=nl(indROM);
        vect=vNL_ROM*(MatInv*nlROM);
    	fROMi(j,k,:)=v_ROM'*(mass*vect);
    end
end
%
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
%in=7500
in=0
THor=7501
fROM=fROMi;

rangeDNSbrut=(in+1):(in+THor);
t0=DNS_t2(rangeDNSbrut);
model_t = t0;

z = zeros(nROM, length(model_t)); % The solution array.
z(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition.

%Use 1st-order accurate time scheme for the first step.
for j = 2
z(:, 2) = (eye(nROM) - dt * aROM) \ (z(:, 1) + dt * dydt_nonlinear(z(:, 1)));
end

%Use 2nd-order accurate time scheme for all subsequent steps.
for j = 3:length(model_t)
    z(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * z(:, j-1) - z(:, j-2) ...
        + 4 * dt * dydt_nonlinear(z(:, j-1))- 2 * dt * dydt_nonlinear(z(:, j-2)));
end

figure(2)
plot(DNS_t,yy(1,:),'r')
hold on
plot(model_t,z(1,:),'k')
hold off

RangeMeas=16;
% 5 6 = (1,0
% 15 16 = (3,0
% 25 26 = (5,0
% 35 36 = (7,0
% 45 46 = (9,0
% 55 56 = (11,0
% 65 66 = (13,0
% 75 76 = (15,0
% 85 86 = (17,0
% 95 96 = (19,0
mDNSbrut=mesDNSbrut(RangeMeas,rangeDNSbrut);
CROM=C(:,1:nROM);
mesROM=CROM*z;
mMODEL=mesROM(RangeMeas,:);
% Compute 2-norm error (on all measurements) on full sequence
errorMeasure=norm(mDNSbrut-mMODEL,'fro')/norm(mDNSbrut,'fro')*100;
%
figure(12)
plot(t0,mDNSbrut(1,:),'r')
hold on
plot(t0,mMODEL(1,:),'b')
hold off
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
