function [c,ceq] = nonlineardynamics(P,D,N,ic0,icf,param)

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

% Get Earth rotation in NED
omega = param.omega;
p = omega(1);
q = omega(2);
r = omega(3);


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

% Get gravity along the trajectory
[gx,gy,gz] = gravityfn(x,y,z,N,param);

% Get aerodynamic drag

rho0 = param.rho0;
h0   = param.h0;
Cd   = param.Cd;
d    = param.d;

rho = rho0*exp(z/h0);
Vsq   = vx.^2 + vy.^2 + vz.^2;
V     = sqrt(Vsq);
FD  = 0.5*rho.*Vsq*Cd*pi*d^2/4;

Fdx = -FD.*(vx./V);
Fdy = -FD.*(vy./V);
Fdz = -FD.*(vz./V);


% Formulate the thrust force
% Main assumption thrust is acting through the body x-axis.

Isp = param.Isp;
g0  = param.g0;

T = 2/tf*D*M*(-Isp*g0);

Tx =  T.*cos(theta).*cos(psi);
Ty =  T.*cos(theta).*sin(psi);
Tz = -T.*sin(theta);


% Equality path constraints
ceq(1:N+1,1)     = 2/tf*D*x - vx;
ceq(N+2:2*N+2)   = 2/tf*D*y - vy;
ceq(2*N+3:3*N+3) = 2/tf*D*z - vz;

ceq(3*N+4:4*N+4) = M.*(2/tf*D*vx + 2*(q*vz - r*vy) + (p*q*y - q^2*x - r^2*x + r*p*z) - gx) - Fdx - Tx;
ceq(4*N+5:5*N+5) = M.*(2/tf*D*vy + 2*(r*vx - p*vz) + (p*q*x + r*q*z - p^2*y - r^2*y) - gy) - Fdy - Ty;
ceq(5*N+6:6*N+6) = M.*(2/tf*D*vz + 2*(p*vy - q*vx) + (p*r*x + p^2*z - q^2*z - r*q*y) - gz) - Fdz - Tz;

ceq(6*N+7:7*N+7) = 2/tf*D*theta - utheta;
ceq(7*N+8:8*N+8) = 2/tf*D*psi   - upsi;

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
bf(8,1) = psi(N+1)   - psif;


ceq = [ceq;bi;bf];

Tmax = param.Tmax;

c(1:N+1,1)     = -T + 1e-10;
c(N+2:2*N+2,1) =  T - Tmax;
% c(2*N+3,1)     = -M(N+1) + mf;

end