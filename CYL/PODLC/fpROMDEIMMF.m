clear;
%%
load('PODMF.mat');
load('PODNLMF.mat');
load('ROMDEIMMF.mat');
load('MatrixFF.mat');
load('coord.mat');
%%
nROM=16;
v_ROM=v(:,1:nROM);
osROM=os(1:nROM);
aROM=a(1:nROM,1:nROM);
eig(aROM)
%%
nROM2=10;
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
resDEIM = @(z) osROM + aROM*z + bROM * ((uuROM * z) .* (qdvdxROM * z)+(vvROM * z) .* (qdvdyROM * z)); 
%%
y = zeros(nROM,1); % The solution array.
%%
deltay=1;
%while norm(deltay)>1e-10
    
    residual=resDEIM(y);
    norm(residual)

    J=zeros(size(aROM));
    for i=1:size(y,1)
        J(:,i)=bROM * (uuROM(:,i) .* (qdvdxROM * y)+(uuROM * y) .* (qdvdxROM(:,i))+vvROM(:,i) .* (qdvdyROM * y)+(vvROM * y) .* (qdvdyROM(:,i)));
    end;
    J=J+aROM;

    deltay=-J\residual;
    norm(deltay)
    y=y+deltay;
%end
%% Load the data.
%%
fid=fopen(strcat('mfcorr.txt'),'rt');
num=fscanf(fid,'%d',1);
b=textscan(fid,'%f',num);
mfcorr=b{1};
fclose(fid);
%%
yy = zeros(nROM,1); % The solution array.
yy=-v_ROM' * mass * mfcorr;
diff=mfcorr+v_ROM*yy;
sqrt(diff'*mass*diff)/sqrt(mfcorr'*mass*mfcorr) % potentiel de reconstruction, doit etre le plus petit possible
%%
diff2=mfcorr+v_ROM*y;
sqrt(diff2'*mass*diff2)/sqrt(mfcorr'*mass*mfcorr)   % reconstruction effective
%%
eig(J)
