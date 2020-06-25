clear
load('Snap_State_POD_MF.mat');
load('ROMGALMF.mat');
%%
% 2nd order time integration, explicit in nonlinear terms, implicit in linear
% terms.
in=15000
model_t = (in+1)*dt:dt:DNS_t(end);
model_nt = length(model_t);
%%
nROM=30;
aROM=a(1:nROM,1:nROM);
osROM=os(1:nROM);
global fROM;
fROM=f(1:nROM,1:nROM,1:nROM);
%%
y = zeros(nROM, model_nt); % The solution array.
y(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition. Try 0 / 1500 / 3000

%Use 1st-order accurate time scheme for the first step.
for j = 2
y(:, 2) = (eye(nROM) - dt * aROM) \ ...y(:,1)
    (y(:, 1) + dt*dydt_nonlinear(y(:, 1)) + dt*osROM);
end

%Use 2nd-order accurate time scheme for all subsequent steps.
for j = 3:model_nt
    y(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * y(:, j-1) - y(:, j-2) ...
        + 4 * dt * dydt_nonlinear(y(:, j-1))- 2 * dt * dydt_nonlinear(y(:, j-2))+ 2 * dt * osROM);
end
%%
figure(2)
plot(DNS_t,yy(1,:),'r')
hold on
plot(model_t,y(1,:),'k')
hold off
%% Error on full snapshot
in=15000
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
fid=fopen('errorvstime_galmf.txt','w')
for i=1:size(timeRange,2)
 fprintf(fid,'%16.8e %16.8e %16.8e\n',timeRange(i),errorMODEL(i),errorTRUNC(i));
end
fclose(fid)
% read measurementmatrix 
%% Read measurements of POD
fid=fopen('mespodmf.txt','rt');
CC=textscan(fid,'%f');
C=reshape(CC{1},100,size(a,1));
clear CC;
fclose(fid);%%
%% Read measurements of mean-flow
fid=fopen('mesmf.txt','rt');
CC=textscan(fid,'%f');
mesMF=CC{1};
clear CC;
fclose(fid);
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
in=15000
rangeDNSbrut=(in+1):step:size(DNS_t2,2);
t0=DNS_t2(rangeDNSbrut);
mDNSbrut=mesDNSbrut(:,rangeDNSbrut);
%% Select same range in DNS measurements from projection
rangeDNS=in/step+1:size(DNS_t,2);
t1=DNS_t(rangeDNS);
mDNS=mesMF*ones(1,size(t1,2))+mesDNS(:,rangeDNS);
%% Check that DNS measurements coincide
figure(12)
numPOS=6;
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
%% Same as above but just for one measurement numPOS
numPOS=5;
CDNS=C(numPOS,1:end);
mesDNS=CDNS*yy;
%%
in=15000
rangeDNSbrut=(in+1):step:size(DNS_t2,2);
t0=DNS_t2(rangeDNSbrut);
mDNSbrut=mesDNSbrut(numPOS,rangeDNSbrut);

rangeDNS=in/step+1:size(DNS_t,2);
t1=DNS_t(rangeDNS);
mDNS=mesMF(numPOS)+mesDNS(rangeDNS);
%%
error=norm(mDNS-mDNSbrut)/norm(mDNSbrut)
%%
figure(12)
plot(t0,mDNSbrut,'r')
hold on
plot(t1,mDNS,'b')
hold off
%%
CROM=C(numPOS,1:nROM);
mesROM=CROM*y;

%%
t0=DNS_t(rangeDNS);
mDNS=mesDNS(rangeDNS);
t1=model_t(rangeMODEL);
mMODEL=mesROM(rangeMODEL);
error=norm(mMODEL-mDNS)/norm(mDNS)
%%
figure(3)
plot(DNS_t,mesDNS(:),'r')
hold on
plot(model_t,mesROM(:),'k')
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
