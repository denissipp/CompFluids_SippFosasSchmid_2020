clear;
%%
load('POD.mat');
load('PODNL.mat');
load('ROMDEIM.mat');
load('MatrixFF.mat');
load('coord.mat');
%%
nROM=6;
v_ROM=v(:,1:nROM);
aROM=a(1:nROM,1:nROM);
eig(aROM)
%%
nROM2=11;
%%
vNL_ROM=vNL(:,1:nROM2);
indROM = ind(1:nROM2);
max(xx(indROM))
bROM = (v_ROM' * mass * vNL_ROM) / vNL_ROM(indROM, :); % Pre-multiplier on nonlinear function.
uuROM = uu(indROM, 1:nROM);
vvROM = vv(indROM, 1:nROM);
qdvdxROM = qdvdx(indROM, 1:nROM); % Truncated pre-multiplier on POD coefficients for dw/dx.
qdvdyROM = qdvdy(indROM, 1:nROM); % Truncated pre-multiplier on POD coefficients for dw/dy.
%%
% The time derivative of the POD coefficients y.
dydt_nonlinear_deim = @(z) bROM * ((uuROM * z) .* (qdvdxROM * z)+(vvROM * z) .* (qdvdyROM * z)); 
%%
%%
load('Snap_State_POD.mat');
%%
%%
% 2nd order time integration, explicit in nonlinear terms, implicit in linear
% terms.
in=15000
model_t = (in+1)*dt:dt:5*DNS_t(end);
model_nt = length(model_t);

%%
z = zeros(nROM, model_nt); % The solution array.
z(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition.

%Use 1st-order accurate time scheme for the first step.
for j = 2
z(:, 2) = (eye(nROM) - dt * aROM) \ (z(:, 1) + dt * dydt_nonlinear_deim(z(:, 1)));
end

%Use 2nd-order accurate time scheme for all subsequent steps.
for j = 3:model_nt
    z(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * z(:, j-1) - z(:, j-2) ...
        + 4 * dt * dydt_nonlinear_deim(z(:, j-1))- 2 * dt * dydt_nonlinear_deim(z(:, j-2)));
end
%%
figure(2)
plot(DNS_t,yy(2,:),'r')
hold on
plot(model_t,z(2,:),'k')
hold off
%%axis([0 400 -5 5])
%%
%%
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
mesROM=CROM*z;
figure(3)
plot(DNS_t,mesDNS(:),'r')
hold on
plot(model_t,mesROM(:),'k')
hold off







