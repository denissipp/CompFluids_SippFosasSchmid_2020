clear
load('Snap_State_POD.mat');
load('ROMGAL.mat');
%%
% 2nd order time integration, explicit in nonlinear terms, implicit in linear
% terms.
in=15000
model_t = (in+1)*dt:dt:2*DNS_t(end);
model_nt = length(model_t);
%%
nROM=10;
aROM=a(1:nROM,1:nROM);
global fROM;
fROM=f(1:nROM,1:nROM,1:nROM);
%%
y = zeros(nROM, model_nt); % The solution array.
y(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition.

%Use 1st-order accurate time scheme for the first step.
for j = 2
y(:, 2) = (eye(nROM) - dt * aROM) \ ...y(:,1)
    (y(:, 1) + dt*dydt_nonlinear(y(:, 1)));
end

%Use 2nd-order accurate time scheme for all subsequent steps.
for j = 3:model_nt
    y(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * y(:, j-1) - y(:, j-2) ...
        + 4 * dt * dydt_nonlinear(y(:, j-1))- 2 * dt * dydt_nonlinear(y(:, j-2)));
end
%%
figure(2)
plot(DNS_t,yy(2,:),'r')
hold on
plot(model_t,y(2,:),'k')
hold off
% read measurementmatrix 
fid=fopen('mespod.txt','rt');
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
