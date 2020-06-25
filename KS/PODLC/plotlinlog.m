
function plotlinlog(name,i,j,fig)

A=load(name);

figure(fig)
semilogy(A(:,i),A(:,j))
hold off;
