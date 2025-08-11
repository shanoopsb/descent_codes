clear;clc;close all;

%% Vehicle Parameters

param.Isp = 282;
param.TLO = 5886*1000;
param.d   = 3.66;
param.Cd  = 0.75;
param.h0  = 7500;
param.g0  = 9.80665;
param.rho0 = 1.2256;
param.t0   = 0;
param.t1   = 40.8;
param.m0   = 76501; 
param.mdry = 25600;

param.T1max = 1/3*param.TLO;
param.T3max = 1/9*param.TLO;

%% 

N1 = 9;
N2 = 7;
N3 = 11;

[tau1,D1] = Dlgr(N1);
[tau2,D2] = Dlgr(N2);
[tau3,D3] = Dlgr(N3);


param.N1 = N1;
param.N2 = N2;
param.N3 = N3;

param.tau1 = tau1;
param.tau2 = tau2;
param.tau3 = tau3;

param.D1 = D1;
param.D2 = D2;
param.D3 = D3;


tspan1 = param.t0 + 0.5*(param.t1 - param.t0)*(1+tau1);
param.tspan1 = tspan1;

%%

function [c,ceq] = nonldynamics(P,ic0,icf,param)

% find P1 ,P2 P3
N1 = param.N1; D1 = param.D1;
N2 = param.N2; D2 = param.D2;
N3 = param.N3; D3 = param.D3;


x1  = P1(1:N1+1);
vx1 = P1(N1+2:2*N1+2);
y1  = P1(2*N1+3:3*N1+3);
vy1 = P1(3*N1+4:4*N1+4);
m1  = P1(4*N1+5:5*N1+5);
k1  = P1(5*N1+6:6*N1+6);

x2  = P2(1:N2+1);
vx2 = P2(N2+2:2*N2+2);
y2  = P2(2*N2+3:3*N2+3);
vy2 = P2(3*N2+4:4*N2+4);
t2  = P2(end); 

x3  = P3(1:N3+1);
vx3 = P3(N3+2:2*N3+2);
y3  = P3(2*N3+3:3*N3+3);
vy3 = P3(3*N3+4:4*N3+4);
m3  = P3(4*N3+5:5*N3+5);
k3  = P3(5*N3+6:6*N3+6);
theta3 = P3(6*N3+7:7*N3+7);
t3     = P3(end);

theta1 = pi;

x0 = ic0(1); vx0 = ic0(2);
y0 = ic0(3); vy0 = ic0(4);
m0 = ic0(5);   
theta30 = ic0(6);


xf = icf(1); vxf = icf(2);
yf = icf(3); vyf = icf(4);
mf = icf(5);
theta3f = icf(6);

tspan1 = param.tspan1;
t0     = param.t0;
t1     = param.t1;
T1     = param.T1max;
T3     = param.T3max;

ceqx1  = 2/(t1-t0)*D1*x1 - vx1;
ceqvx1 = m1.*(2/(t1-t0)*D1*vx1) - T1*k1.*cos(theta1) + FD1.*(vx1./V1);
ceqy1  = 2/(t1-t0)*D1*y1 - vy1;
ceqvy1 = m1.*(2/(t1-t0)*D1*vy1 + g0) - T1*k1.*sin(theta1) + FD1.*(vy1./V1);
ceqm1  = (2/(t1-t0)*D1*m1*Isp*g0) + T1*k1;

bi(1,1) = x1(1)  - x0;
bi(2,1) = vx1(1) - vx0;
bi(3,1) = y1(1)  - y0;
bi(4,1) = y1(1)  - vy0;
bi(5,1) = m1(1) - m0;


ceq1 = [ceqx1;ceqvx1;ceqy1;ceqvy1;ceqm1;bi];

Li1(1,1) = x1(N1+1) - x2(1);
Li1(2,1) = vx1(N1+1) - vx2(1);
Li1(3,1) = y1(N1+1) - y2(1);
Li1(4,1) = vy1(N1+1) - vy2(1);

m2 = m1(end);

ceqx2  = 2/(t2-t1)*D2*x2 - vx2;
ceqvx2 = m2.*(2/(t2-t1)*D2*vx2) - FD2.*(vx2./V2);
ceqy2  = 2/(t2-t1)*D2*y2 - vy2;
ceqvy2 = m2.*(2/(t2-t1)*D2*vy2 + g0) - FD2.*(vy2./V2);

ceq2 = [Li1;ceqx2;ceqvx2;ceqy2;ceqvy2];

Li2(1,1) = x2(N2+1) - x3(1);
Li2(2,1) = vx2(N2+1) - vx3(1);
Li2(3,1) = y2(N2+1) - y3(1);
Li2(4,1) = vy2(N2+1) - vy3(1);
Li2(5,1) = m2 - m3(1);
Li2(6,1) = theta3(1) - theta30;


ceqx3  = 2/(t3-t2)*D3*x3 - vx3;
ceqvx3 = m3.*(2/(t3-t2)*D3*vx3) - T3*k3.*cos(theta3) - FD3.*(vx3./V3);
ceqy3  = 2/(t3-t2)*D3*y3 - vy3;
ceqvy3 = m3.*(2/(t3-t2)*D3*vy3 + g0) - T3*k3.*sin(theta3) - FD3.*(vy3./V3);
ceqm3  = (2/(t3-t2)*D3*m3*Isp*g0) + T3*k3;

bf(1,1) = x3(N3+1) - xf;
bf(2,1) = vx3(N3+1) - vxf;
bf(3,1) = y3(N3+1) - yf;
bf(4,1) = vy3(N3+1) - vyf;
bf(6,1) = theta3(N3+1) - theta3f;

ceq3 = [Li2;ceqx3;ceqvx3;ceqy3;ceqvy3;ceqm3;bf];

ceq = [ceq1;ceq2;ceq3];

ck1m = -k1;
ck1p =  k1 - 1;

ck3m = -k3;
ck3p =  k3 - 1;




end


