clear
load('ROMGALMF.mat');
%%
%%
nROM=60;
global aROM;
global osROM;
global fROM;
aROM=a(1:nROM,1:nROM);
osROM=os(1:nROM);
fROM=f(1:nROM,1:nROM,1:nROM);
%%
y = zeros(nROM,1); % The solution array.
%%
deltay=1;
while norm(deltay)>1e-13
    
    residual=res(y);
    norm(residual)

    J=zeros(size(aROM));
    for i=1:size(y,1)
        JJ=fROM(:,:,i);
        J(i,:)=y'*(JJ+JJ');
    end;
    J=J+aROM;

    deltay=-J\residual;
    norm(deltay)
    y=y+deltay;
end
%% Load the data.
load('PODMF.mat');
load('MatrixFF.mat');
v_ROM=v(:,1:nROM);
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
%%
y
