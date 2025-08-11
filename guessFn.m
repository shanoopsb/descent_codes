function Pg = guessFn(D,N,tau,ic0,param)

y0 = ic0(1:7);
% U   = ic0(8:9);

U = [deg2rad(90);deg2rad(135)];

tspan = 0:0.01:300;
odefun = @(t,Y) dynamics_guess(t,Y,U,param);


[tout,yout] = ode45(odefun,tspan,y0);

x = yout(:,1);
y = yout(:,1);
z = yout(:,1);

vx = yout(:,1);
vy = yout(:,1);
vz = yout(:,1);

m = yout(:,1);

Tguess = tout(end)/2*(1+tau);

xg = spline(tout,x,Tguess);
yg = spline(tout,y,Tguess);
zg = spline(tout,z,Tguess);


vxg = spline(tout,vx,Tguess);
vyg = spline(tout,vy,Tguess);
vzg = spline(tout,vz,Tguess);

mg = spline(tout,m,Tguess);

thetag = U(1)*ones(N+1,1);
psig   = U(2)*ones(N+1,1);

uthetag = 2/tout(end)*D*thetag;
upsig   = 2/tout(end)*D*psig;

tfg = tout(end);

Pg = [xg;yg;zg;vxg;vyg;vzg;mg;thetag;psig;uthetag;upsig;tfg];

end

function dy = dynamics_guess(t,Y,U,param)

x  = Y(1);
y  = Y(2);
z  = Y(3);
vx = Y(4);
vy = Y(5);
vz = Y(6);
m  = Y(7);

theta = U(1);
psi   = U(2);

if t <= 200
   T = 0;
else 
    T = param.Tmax;
end

rT    = param.rT;
omega = param.omega;
r = [x;y;z];
V = [vx;vy;vz];

rho0 = param.rho0;
h0   = param.h0;
Cd   = param.Cd;
d    = param.d;

rho      = rho0*exp(z/h0);
Vsq      = vx^2 + vy^2 + vz^2;
% Vmag     = sqrt(Vsq);
FD       = 0.5*rho*Vsq*Cd*pi*d^2/4;

Fdx = -FD*(vx/norm(V));
Fdy = -FD*(vy/norm(V));
Fdz = -FD*(vz/norm(V));

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

dy(4:6,1) =  a_thrust + a_aero + a_grav - 2*cross(omega,V) - cross(omega,cross(omega,r)) ;
dy(7,1)   = -T/(param.Isp*param.g0);
end


