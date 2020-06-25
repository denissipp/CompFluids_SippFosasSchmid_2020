clear
load('Snap_State_POD_MF.mat');
load('ROMGALMF.mat');
%%
load('mes_0.mat')
%% Read measurements of POD
fid=fopen('mespodmf.txt','rt');
CC=textscan(fid,'%f');
C=reshape(CC{1},61,size(a,1));
clear CC;
fclose(fid);%%
%% Read measurements of mean-flow
fid=fopen('mesmf.txt','rt');
CC=textscan(fid,'%f');
mesMF=CC{1};
clear CC;
fclose(fid);
%%
in=20000
step=20;
rangeDNSbrut=(in+1):size(DNS_t2,2);
t0=DNS_t2(rangeDNSbrut);

model_t = (in+1)*dt:dt:DNS_t2(end);
rangeMODEL=1:(DNS_t2(end)-model_t(1))/dt+1;
t1=model_t(rangeMODEL);
%%
nROM=8;
aROM=a(1:nROM,1:nROM);
osROM=os(1:nROM);
global fROM;
fROM=f(1:nROM,1:nROM,1:nROM);
%%
y = zeros(nROM, length(model_t)); % The solution array.
y(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition. Try 0 / 1500 / 3000
%%
%Use 1st-order accurate time scheme for the first step.
for j = 2
y(:, 2) = (eye(nROM) - dt * aROM) \ ...y(:,1)
    (y(:, 1) + dt*dydt_nonlinear(y(:, 1)) + dt*osROM);
end

%Use 2nd-order accurate time scheme for all subsequent steps.
for j = 3:length(model_t)
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
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%% ERRORS ON ROM MEAS
%% Select a range, and swith to Step=5 for DNS measurements
%RangeMeas=1:2:100;
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
mesROM=CROM*y;
mMODEL=mesROM(RangeMeas,rangeMODEL)+mesMF(RangeMeas)*ones(1,size(rangeMODEL,2));
%%
figure(12)
plot(t0,mDNSbrut(1,:),'r')
hold on
plot(t1,mMODEL(1,:),'b')
hold off
%%
%%
fid=fopen('meascomp_mf_25.txt','w')
for i=1:size(t0,2)
 fprintf(fid,'%16.8e %16.8e\n',t0(i),mDNSbrut(i));
end
fclose(fid)
%%
fid=fopen('meascomp_gal_mf_25.txt','w')
for i=1:size(t0,2)
 fprintf(fid,'%16.8e %16.8e\n',t1(i),mMODEL(i));
end
fclose(fid)



%%



% read measurementmatrix 
fid=fopen('mespodmf.txt','rt');
CC=textscan(fid,'%f');
C=reshape(CC{1},100,size(a,1));
clear CC;
fclose(fid);
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

%%
numPOS=95;
CDNS=C(numPOS,1:end);
mesDNS=CDNS*yy;
CROM=C(numPOS,1:nROM);
mesROM=CROM*y;
figure(3)
plot(DNS_t,mesDNS(:),'r')
hold on
plot(model_t,mesROM(:),'k')
hold off
