clear;
%%
%%
load('MatrixFF.mat');
load('SnapBF.mat');

%%
METHOD='DEIM';
TRAJ='BF';
disp(['METHOD=',METHOD,' TRAJ=',TRAJ])
nROM=30;
nROM2=20;

file = strcat(METHOD,'-',num2str(nROM),'-',num2str(nROM2),'-',TRAJ,'.txt');
fid = fopen(file,'w');

disp(['nROM=', num2str(nROM), ' nROM2=',num2str(nROM2)])

if TRAJ=='BF'
    load('POD.mat');
    load('PODNL.mat');
    load('ROMDEIM.mat');
    load('Snap_State_POD.mat');
else
    load('PODMF.mat');
    load('PODNLMF.mat');
    load('ROMDEIMMF.mat');
    load('Snap_State_POD_MF.mat');
end

%%
v_ROM=v(:,1:nROM);
if TRAJ=='MF'
	osROM=os(1:nROM);
else
	osROM=zeros(nROM,1);
end

aROM=a(1:nROM,1:nROM);
[evectL,ev]=eig(aROM);
maxvp=-1e30;
for ii=1:nROM
    if (real(ev(ii,ii))>maxvp)&&(imag(ev(ii,ii))>=0)
        maxvp=real(ev(ii,ii));
        indeL=ii;
    end
end
vect1= v_ROM*evectL(:,indeL);
if TRAJ=='BF'
	vect2 =eigen;
	tev = 0.1307763186079424+0.8173834768895367*1j;
else
	vect2 =eigenMF;
	tev= 0.001949689325986477+1j*1.059033061124873;
end
psLIN = abs(vect1'*(mass*vect2))/sqrt(abs(vect1'*(mass*vect1)))/sqrt(abs(vect2'*(mass*vect2)));
erreigLIN=abs(ev(indeL,indeL)-tev)/abs(tev)*100;
disp(['erreigLIN=', num2str(erreigLIN,2), ' psLIN=', num2str((1-psLIN)*100,2)])
fprintf(fid,'erreigLIN=%5.2f / psLIN=%5.2f\n\n',erreigLIN,(1-psLIN)*100);
%%
fROMi=zeros(nROM,nROM,nROM);
%% GAL1
if METHOD=='GAL1'
    for j=1:nROM
        for k=1:nROM
            vect=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j)+vv(:,j).*qdvdy(:,k)+vv(:,k).*qdvdy(:,j));
            fROMi(j,k,:)=v_ROM'*(mass*vect);   
        end
    end
elseif METHOD=='GAL2'
        %% GAL2
        vNL_ROM=vNL(:,1:nROM2);
        for j=1:nROM
            for k=1:nROM
                nl=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j)+vv(:,j).*qdvdy(:,k)+vv(:,k).*qdvdy(:,j));
                vect=vNL_ROM*(vNL_ROM'*(mass*nl));
                fROMi(j,k,:)=v_ROM'*(mass*vect);
            end
        end
else
    %%DEIM
    vNL_ROM=vNL(:,1:nROM2);
    indROM = ind(1:nROM2);
    MatInv=pinv(vNL_ROM(indROM, :));
    %%
    fROMi=zeros(nROM,nROM,nROM);
    for j=1:nROM
        for k=1:nROM
            nl=0.5*(uu(:,j).*qdvdx(:,k)+uu(:,k).*qdvdx(:,j)+vv(:,j).*qdvdy(:,k)+vv(:,k).*qdvdy(:,j));
            nlROM=nl(indROM);
            vect=vNL_ROM*(MatInv*nlROM);
            fROMi(j,k,:)=v_ROM'*(mass*vect);
        end
    end
end

%%
global fROM;
fROM=fROMi;
%%
if TRAJ=='MF'
	x0=zeros(nROM,1);
	for ii=1:20
    		res=osROM+aROM*x0+dydt_nonlinear(x0);
    		mat=aROM;
    		for i=1:nROM
        		for j=1:nROM
            			for k=1:nROM
                			mat(i,j)=mat(i,j)+(fROM(j,k,i)*x0(k)+fROM(k,j,i)*x0(k));
            			end
        		end
    		end
    		dx0=-mat\res;
    		x0=x0+dx0;
	end
	diff = v(:,1:nROM)*x0-mfcorr;
	errbf=sqrt((diff'*(mass*diff))/(mfcorr'*(mass*mfcorr)))*100.;
	disp(['errbf=', num2str(errbf,2)])
    fprintf(fid,'errbf=%5.2f\n',errbf);

	[evect,ev]=eig(mat);
	maxvp=-1e30;
	nu=0;
	for ii=1:nROM
    		if real(ev(ii,ii))>0
        		nu=nu+1;
    		end
    		if (real(ev(ii,ii))>maxvp)&&(imag(ev(ii,ii))>=0)
        		maxvp=real(ev(ii,ii));
        		inde=ii;
    		end
	end
	vect1 = v_ROM*evect(:,inde);
	vect2 = eigen;
	ps = abs(vect1'*(mass*vect2))/sqrt(abs(vect1'*(mass*vect1)))/sqrt(abs(vect2'*(mass*vect2)));
	tev = 0.1307763186079424+0.8173834768895367*1j;
	erreigbf=abs(ev(inde,inde)-tev)/abs(tev)*100;
	disp(['nuBF=', num2str(nu), ' erreigBF=',num2str(erreigbf,2), ' psBF=', num2str((1-ps)*100,2)])
    fprintf(fid,'nuBF=%d / erreigBF=%5.2f / psBF=%5.2f\n\n',nu,erreigbf,(1-ps)*100);
end
    %%
    meas = zeros(15000,1)
    z = zeros(nROM, 1); % The solution array.
    if TRAJ=='BF'
        z(1:nROM,1)=-v_ROM'*mass*mfcorr;
    else
        z(1:nROM,1)=v_ROM'*mass*mfcorr;
    end
    zp = z;
    z(:) = (eye(nROM) - dt * aROM) \ (zp(:) + dt * dydt_nonlinear(zp(:)) + dt * osROM);   
    for j = 1:15000
        zpp=zp;
        zp=z;
        z(:) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
            (4 * zp(:) - zpp(:) ...
            + 4 * dt * dydt_nonlinear(zp(:))- 2 * dt * dydt_nonlinear(zpp(:)) + 2 * dt * osROM);
        meas(j)=z(1);
    end
    %%
    zm = z;
    for j = 1:15000
        zpp=zp;
        zp=z;
        z(:) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
            (4 * zp(:) - zpp(:) ...
            + 4 * dt * dydt_nonlinear(zp(:))- 2 * dt * dydt_nonlinear(zpp(:)) + 2 * dt * osROM);
        zm = (j*zm+z)/(j+1);
    end
    
   %%
   if TRAJ=='BF'
       diff = v(:,1:nROM)*zm+mfcorr;
   else
       diff = v(:,1:nROM)*zm-mfcorr;
   end
    errmf=sqrt((diff'*(mass*diff))/(mfcorr'*(mass*mfcorr)))*100.;
	disp(['errmf=', num2str(errmf,2)])
    fprintf(fid,'errmf=%5.2f\n',errmf);

%%
    mat=aROM;
    for i=1:nROM
        for j=1:nROM
   			for k=1:nROM
               			mat(i,j)=mat(i,j)+(fROM(j,k,i)*zm(k)+fROM(k,j,i)*zm(k));
            end
        end
    end

	[evect,ev]=eig(mat);
	maxvp=-1e30;
	nu=0;
	for ii=1:nROM
    		if real(ev(ii,ii))>0
        		nu=nu+1;
    		end
    		if (real(ev(ii,ii))>maxvp)&&(imag(ev(ii,ii))>=0)
        		maxvp=real(ev(ii,ii));
        		inde=ii;
    		end
	end
	vect1 = v_ROM*evect(:,inde);
	vect2 = eigenMF;
	ps = abs(vect1'*(mass*vect2))/sqrt(abs(vect1'*(mass*vect1)))/sqrt(abs(vect2'*(mass*vect2)));
	tev = 0.001949689325986477+1j*1.059033061124873;
	erreigMF=abs(ev(inde,inde)-tev)/abs(tev)*100;
	disp(['nuMF=', num2str(nu), ' erreigMF=',num2str(erreigMF,2), ' psMF=', num2str((1-ps)*100,2)])
    fprintf(fid,'nuMF=%d / erreigMF=%5.2f / psMF=%5.2f\n\n',nu,erreigMF,(1-ps)*100);


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
assymetry=sqrt(num/den)*100;

disp(['assym=',num2str(assymetry,2)])

%%

for kkk=0:1
    if kkk==0
            in=0; disp 'TR';    fprintf(fid,'TR\n');
    else
            in=7500; disp 'LC'; fprintf(fid,'\nLC\n');
    end
        
    for lll=0:1
        if lll==0
            fROM=fROMi; disp 'fROMi'; fprintf(fid,'fROMi : '); 
        else
            fROM=fROMa; disp 'fROMa'; fprintf(fid,'fROMa : ');
        end

        
THor=7501;

rangeDNSbrut=(in+1):(in+THor);
t0=rangeDNSbrut*dt;
model_t = t0;

z = zeros(nROM, length(model_t)); % The solution array.
z(1:nROM, 1) = yy(1:nROM,in/step+1); % The inital condition.
z(:, 2) = (eye(nROM) - dt * aROM) \ (z(:, 1) + dt * dydt_nonlinear(z(:, 1)) + dt * osROM);
for j = 3:length(model_t)
    z(:, j) = (3 * eye(nROM) - 2 * dt * aROM) \ ...
        (4 * z(:, j-1) - z(:, j-2) ...
        + 4 * dt * dydt_nonlinear(z(:, j-1))- 2 * dt * dydt_nonlinear(z(:, j-2)) + 2 * dt * osROM);
end
rangeDNS=in/step+1:(in/step+THor/step+1);
rangeMODEL=1:step:THor;
yye=yy(:,rangeDNS);
ye=z(:,rangeMODEL);
A=norm(ye(1:nROM,:)-yye(1:nROM,:),'fro')^2/size(yye,2); % mean error per snapshot
B=norm(yye(nROM+1:end,:),'fro')^2/size(yye,2);
CC=norm(yye(:,:),'fro')^2/size(yye,2);
errMODEL=sqrt(A/CC)*100;
errTRUNC=sqrt(B/CC)*100;
    
disp([num2str(errTRUNC,2), ' ', num2str(errMODEL,2)]);
fprintf(fid,'errT=%5.2f / errM=%5.2f\n',errTRUNC,errMODEL);

    end
end
fclose(fid);
