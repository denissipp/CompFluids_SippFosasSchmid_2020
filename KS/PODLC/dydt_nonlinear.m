function rhs=dydt_nonlinear(y)
global fROM;
rhs=zeros(size(y));
for i=1:size(y,1)
    rhs(i)=y'*fROM(:,:,i)*y;
end;
end
