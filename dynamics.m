function dy = dynamics(Y,U,param)

x  = Y(1);
y  = Y(2);
z  = Y(3);
vx = Y(4);
vy = Y(5);
vz = Y(6);
m  = Y(7);

T     = U(1);
theta = U(2);
psi   = U(3);

rT    = param.rT;
omega = param.omega;
r = [x;y;z];
V = [vx;vy;vz];

rho0 = param.rho0;
h0   = param.h0;
Cd   = param.Cd;
d    = param.d;

rho      = rho0*exp(z/h0);
Vsq      = vx.^2 + vy.^2 + vz.^2;
% Vmag     = sqrt(Vsq);
FD       = 0.5*rho.*Vsq*Cd*pi*d^2/4;

Fdx = -FD.*(vx./norm(V));
Fdy = -FD.*(vy./norm(V));
Fdz = -FD.*(vz./norm(V));

a_aero   = 1/m*[Fdx;Fdy;Fdz];

Mu = param.Mu;
a_grav   = -Mu*(rT+r)/(norm(r+rT))^3;

Tx =  T*cos(theta)*cos(psi);
Ty =  T*cos(theta)*sin(psi);
Tz = -T*sin(theta);

a_thrust = 1/m*[Tx;Ty;Tz];

dy(1,1) = vx;
dy(2,1) = vy;
dy(3,1) = vz;

dy(4:6,1) = a_thrust + a_aero + a_grav - 2*cross(omega,V) - cross(omega,cross(omega,r)) ;
dy(7,1)   = -T/(param.Isp*param.g0);
end