
function plotlinlin(name,i,j,fig)

A=load(name);

figure(fig)
plot(A(:,i),A(:,j),'x')
hold off;
