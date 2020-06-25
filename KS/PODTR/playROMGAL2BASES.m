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
in=0
step=1;
rangeDNSbrut=(in+1):step:size(DNS_t2,2);
t0=DNS_t2(rangeDNSbrut);

model_t = (in+1)*dt:dt:DNS_t2(end);
rangeMODEL=1:step:(DNS_t2(end)-model_t(1))/dt+1;
t1=model_t(rangeMODEL);
%%
nROM=40;
v_ROM=v(:,1:nROM);
aROM=a(1:nROM,1:nROM);
eig(aROM)
%%
nROM2=20; % 10/11 OK
%%
vNL_ROM=vNL(:,1:nROM2);
%%
global fROM;
fROM=zeros(nROM,nROM,nROM);
for j=1:nROM
	j
    for k=1:nROM
		nl=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j));
%		nl=uu(:,j).*qdvdx(:,k);
        vect=vNL_ROM*(vNL_ROM'*(mass*nl));
%		for i=1:nROM
%			fROM(j,k,i)=v_ROM(:,i)'*vect(:);
%        end
    	fROM(j,k,:)=v_ROM'*(mass*vect);
    end
end
%%
fROMs=zeros(nROM,nROM,nROM);
for j=1:nROM
	j
    for k=1:nROM
		for i=1:nROM
			fROMs(j,k,i)=(fROM(i,j,k)+fROM(i,k,j)+fROM(j,i,k)+fROM(j,k,i)+fROM(k,i,j)+fROM(k,j,i))/6;
        end
    end
end
%%
fROMa=zeros(nROM,nROM,nROM);
for j=1:nROM
	j
    for k=1:nROM
		for i=1:nROM
			fROMa(j,k,i)=(5*fROM(j,k,i)-fROM(j,i,k)-fROM(i,j,k)-fROM(i,k,j)-fROM(k,i,j)-fROM(k,j,i))/6;
        end
    end
end

%%
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
%%
figure(2)
plot(DNS_t,yy(1,:),'r')
hold on
plot(model_t,z(1,:),'k')
hold off

%%
%%%%%%%%%%%%%%%%%%%%%%%%% ERRORS ON ROM MEAS
%% Select a range, and swith to Step=5 for DNS measurements
%RangeMeas=1:2:61;
RangeMeas=11;
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
%%
mDNSbrut=mesDNSbrut(RangeMeas,rangeDNSbrut);
%%
CROM=C(:,1:nROM);
mesROM=CROM*z;
mMODEL=mesROM(RangeMeas,rangeMODEL);
%%
figure(12)
plot(t0,mDNSbrut(1,:),'r')
hold on
plot(t1,mMODEL(1,:),'b')
hold off
%%
%%
fid=fopen('meascomp_bf_25.txt','w')
for i=1:size(t0,2)
 fprintf(fid,'%16.8e %16.8e\n',t0(i),mDNSbrut(i));
end
fclose(fid)
%%
fid=fopen('meascomp_gal_bf_25.txt','w')
for i=1:size(t0,2)
 fprintf(fid,'%16.8e %16.8e\n',t1(i),mMODEL(i));
end
fclose(fid)

%% Compute 2-norm error (on all measurements) on full sequence
error2=norm(mDNSbrut-mMODEL,'fro')/norm(mDNSbrut,'fro')
%% Analyse 2-norm error (on all measurements) as a function of time
error=zeros(1,size(mDNSbrut,2));
normalfac=norm(mDNSbrut,'fro')/sqrt(size(mDNSbrut,2)*size(mDNSbrut,1));
for i=1:size(mDNSbrut,2)
    error(i)=norm(mMODEL(:,i)-mDNSbrut(:,i))/sqrt(size(mDNSbrut,1))/normalfac;
end
figure(13)
semilogy(error)
%%
fid=fopen('errormeasvstime_gal.txt','w')
for i=1:size(t0,2)
 fprintf(fid,'%16.8e %16.8e\n',t0(i),error(i));
end
fclose(fid)
%% Error on full snapshot
in=0
rangeDNS=in/step+1:size(DNS_t,2);
rangeMODEL=1:step:(DNS_t(end)-model_t(1))/dt+1;
yye=yy(:,rangeDNS);
ye=y(:,rangeMODEL);
timeRange=model_t(rangeMODEL);
%%
nT=nROM; % or 6
A=norm(ye(1:nT,:)-yye(1:nT,:),'fro')^2/size(yye,2); % mean error per snapshot
B=norm(yye(nT+1:end,:),'fro')^2/size(yye,2);
CC=norm(yye(:,:),'fro')^2/size(yye,2);
errMODEL=A/CC
errTRUNC=B/CC
errorMODEL=zeros(1,size(yye,2));
errorTRUNC=zeros(1,size(yye,2));
for i=1:size(yye,2)
    errorMODEL(i)=(norm(ye(1:nT,i)-yye(1:nT,i)))/sqrt(CC);
    errorTRUNC(i)=norm(yye(nT+1:end,i))/sqrt(CC);
end
errMODELbis=norm(errorMODEL)^2/size(yye,2) % recover global mean error
%%
figure(14)
semilogy(timeRange,errorMODEL,'r')
hold on;
semilogy(timeRange,errorTRUNC,'b')
hold off;
%%
fid=fopen('errorvstime_gal.txt','w')
for i=1:size(timeRange,2)
 fprintf(fid,'%16.8e %16.8e %16.8e\n',timeRange(i),errorMODEL(i),errorTRUNC(i));
end
fclose(fid)
%%
%% Read measurements of POD
fid=fopen('mespod.txt','rt');
CC=textscan(fid,'%f');
C=reshape(CC{1},100,size(a,1));
clear CC;
fclose(fid);%%
%% Read measurements of DNS Step=1
fid=fopen('../DNS2/mes_0.txt','rt');
CC=textscan(fid,'%f');
M=reshape(CC{1},101,size(CC{1},1)/101);
DNS_t2=M(1,:);
mesDNSbrut=M(2:101,:);
clear CC,M;
fclose(fid);%%
%%
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
%% Generate DNS measurements from projections on snapshot Step=5
CDNS=C(:,1:end);
mesDNS=CDNS*yy;
%% Select a range, and swith to Step=5 for DNS measurements
in=0
rangeDNSbrut=(in+1):step:size(DNS_t2,2);
t0=DNS_t2(rangeDNSbrut);
mDNSbrut=mesDNSbrut(:,rangeDNSbrut);
%% Select same range in DNS measurements from projection
rangeDNS=in/step+1:size(DNS_t,2);
t1=DNS_t(rangeDNS);
mDNS=mesDNS(:,rangeDNS);
%% Check that DNS measurements coincide
figure(12)
numPOS=5;
plot(t0,mDNSbrut(numPOS,:),'r')
hold on
plot(t1,mDNS(numPOS,:),'b')
hold off
%% Compute 2-norm error (on all measurements) on full sequence
error2=norm(mDNS-mDNSbrut,'fro')/norm(mDNSbrut,'fro')
%% Analyse 2-norm error (on all measurements) as a function of time
error=zeros(1,size(mDNS,2));
for i=1:size(mDNS,2)
    error(i)=norm(mDNS(:,i)-mDNSbrut(:,i))/sqrt(size(mDNSbrut,1));
end
errror=error/(norm(mDNSbrut,'fro')/sqrt(size(mDNSbrut,2)*size(mDNSbrut,1)));
figure(13)
semilogy(error)
%%
figure(12)
plot(t0,mDNSbrut,'r')
hold on
plot(t1,mDNS,'b')
hold off

%%
fid=fopen('compmes.txt','w')
%%
for i=1:size(DNS_t,2)
 fprintf(fid,'%16.8e %16.8e\n',DNS_t(i),mesDNS(i));
end
fclose(fid)
%%
fid=fopen('compmes2.txt','wt');
for i=1:size(model_t,2)
 fprintf(fid,'%16.8e %16.8e\n',model_t(i),mesROM(i));
end
%%
