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
a_grav = gravityfn(x,y,z,N,param);

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




end