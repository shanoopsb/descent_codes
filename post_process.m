function post_process(Popt,N,D,tau,param)


% Extract states
x  = Popt(1:N+1);
y  = Popt(N+2:2*N+2);
z  = Popt(2*N+3:3*N+3);
vx = Popt(3*N+4:4*N+4);
vy = Popt(4*N+5:5*N+5);
vz = Popt(5*N+6:6*N+6);

M = Popt(6*N+7:7*N+7);

theta = Popt(7*N+8:8*N+8);
psi   = Popt(8*N+9:9*N+9);

% Extract control
utheta = Popt(9*N+10:10*N+10);
upsi   = Popt(10*N+11:11*N+11);

tf     = Popt(end);

t0   = param.t0;
Tout = t0 + (tf-t0)/2*(1+tau);

Isp = param.Isp;
g0  = param.g0;

T = 2/tf*D*M*(-Isp*g0);


figure(1)
plot(Tout,x/1000,'r-.','DisplayName','x')
hold on
plot(Tout,y/1000,'g-.','DisplayName','y')
plot(Tout,abs(z)/1000,'b-.','DisplayName','abs(z)')
xlabel("Time (s)")
ylabel("position (km)")
legend("Location","best")

figure(2)
plot(Tout,vx,'r-.','DisplayName','vx')
hold on
plot(Tout,vy,'g-.','DisplayName','vy')
plot(Tout,vz,'b-.','DisplayName','vzv')
xlabel("Time (s)")
ylabel("velocity (m/s)")
legend("Location","best")

figure(3)
plot(Tout,rad2deg(theta),'r-.','DisplayName','\theta')
hold on
plot(Tout,rad2deg(psi),'b-.','DisplayName','\psi')
xlabel("Time (s)")
ylabel("\theta & \psi (deg)")
legend("Location","best")

figure(4)
plot(Tout,rad2deg(utheta),'r-.','DisplayName','u_{\theta}')
hold on
plot(Tout,rad2deg(upsi),'b-.','DisplayName','u_{\psi}')
xlabel("Time (s)")
ylabel("u_{\theta} & u_{\psi} (deg/s)")
legend("Location","best")

figure(5)
plot(Tout,M,'r-.','DisplayName',"Mass")
xlabel("Time (s)")
ylabel("Mass (kg)")
legend("Location","best")

figure(6)
plot(Tout,T/1000,'r-.','DisplayName',"Thrust")
xlabel("Time (s)")
ylabel("Thrust (kN)")
legend("Location","best")



end