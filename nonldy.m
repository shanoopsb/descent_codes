function [c,ceq] = nonldy(P,D,N,ic0,icf,param)

% Extract states
x  = P(1:N+1);
y  = P(N+2:2*N+2);
z  = P(2*N+3:3*N+3);
vx = P(3*N+4:4*N+4);
vy = P(4*N+5:5*N+5);
vz = P(5*N+6:6*N+6);

M = P(6*N+7:7*N+7);

theta = P(7*N+8:8*N+8);
psi   = P(8*N+9:9*N+9);

% Extract control
utheta = P(9*N+10:10*N+10);
upsi   = P(10*N+11:11*N+11);

tf = P(end);

Isp = param.Isp;
g0  = param.g0;
T = 2/tf*D*M*(-Isp*g0);

for i = 1:N+1
    Y = [x(i);y(i);z(i);vx(i);vy(i);vz(i);M(i)];
    U = [T(i);theta(i);psi(i)];

    dy = dynamics(Y,U,param);

    xdot(i,1)  = dy(1);
    ydot(i,1)  = dy(2);
    zdot(i,1)  = dy(3);
    vxdot(i,1) = dy(4);
    vydot(i,1) = dy(5);
    vzdot(i,1) = dy(6);
    mdot(i,1) =  dy(7);
end

ceqx  = 2/tf*D*x - xdot;
ceqy  = 2/tf*D*y - ydot;
ceqz  = 2/tf*D*z - zdot;
ceqvx = 2/tf*D*vx - vxdot;
ceqvy = 2/tf*D*vy - vydot;
ceqvz = 2/tf*D*vz - vzdot;

ceqm  = 2/tf*D*M - mdot;

cequt = 2/tf*D*theta - utheta;
cequp = 2/tf*D*psi   - upsi;

ceq = [ceqx;ceqy;ceqz;ceqvx;ceqvy;ceqvz;ceqm;cequt;cequp];

% Extract initial conditions
x0   = ic0(1); y0   = ic0(2); z0  = ic0(3);
vx0  = ic0(4); vy0  = ic0(5); vz0 = ic0(6);
m0   = ic0(7);
theta0 = ic0(8); psi0 = ic0(9);

% Extract final conditions
xf   = icf(1); yf   = icf(2); zf  = icf(3);
vxf  = icf(4); vyf  = icf(5); vzf = icf(6);
mf   = icf(7);
thetaf = icf(8); psif = icf(9);

% initial conditions
bi(1,1) = x(1)   - x0;
bi(2,1) = y(1)   - y0;
bi(3,1) = z(1)   - z0;
bi(4,1) = vx(1)  - vx0;
bi(5,1) = vy(1)  - vy0;
bi(6,1) = vz(1)  - vz0;
bi(7,1) = M(1)   - m0;
bi(8,1) = theta(1) - theta0;
bi(9,1) = psi(1)   - psi0;

% final conditions
bf(1,1) = x(N+1)   - xf;
bf(2,1) = y(N+1)   - yf;
bf(3,1) = z(N+1)   - zf;
bf(4,1) = vx(N+1)  - vxf;
bf(5,1) = vy(N+1)  - vyf;
bf(6,1) = vz(N+1)  - vzf;
bf(7,1) = theta(N+1) - thetaf;
% bf(8,1) = psi(N+1)   - psif;

ceq = [ceq;bi;bf];


Tmax = param.Tmax;
c(1:N+1,1)     = -T + 1e-10;
c(N+2:2*N+2,1) =  T - Tmax;

end