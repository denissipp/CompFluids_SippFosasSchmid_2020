
function plotlinlin2(name,fig)

A=load(name);

figure(fig)
plot(A(2:end),'x')
hold off;
