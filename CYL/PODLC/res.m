function sol=res(y)
global fROM;
global osROM;
global aROM;
sol=zeros(size(y));
sol=osROM+aROM*y;
for i=1:size(y,1)
    sol(i)=sol(i)+y'*fROM(:,:,i)*y;
end;
end
