


fun = @(t) sin(t);

N = 50;

[tau,D] = Dlgr(N);

t0 = 0;
tf = 2*pi;

tspan = t0 + 0.5*(tf-t0)*(1+tau);
fval = fun(tspan);
fdval = 2/(tf-t0)*D*fval;

f_real = cos(tspan);

plot(tspan,fdval,'ro',DisplayName='LGR')
hold on
plot(tspan,f_real,'b*',DisplayName='Analytical')
xlabel('rad')
ylabel("y")
legend()
