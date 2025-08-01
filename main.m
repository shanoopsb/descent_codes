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

param.Tmax = 1/9*param.TLO;

%% Landing point 
Latitude  = 13.5111;
Longitude = 81.5412;

 
lat  = deg2rad(Latitude);
long = deg2rad(Longitude);

Omega = [0;0;7.29211*10^-5];       % rotation of Earth in ECI

Reci2ned = eci2tcr(lat,long);
omega    = Reci2ned*Omega;         % rotation of Earth in NED

param.omega = omega;

ft2m  = 0.3048;                    % Feet to meters 

RE =  2.0925741 * 10^7 * ft2m;
RP =  2.0855590 * 10^7 * ft2m;
Mu =  1.4076539 * 10^16 * ft2m^3;

K = (RE/RP)^2;
lat_c = atan(tan(lat)/K);
Rs = RE/sqrt(1 + (K-1)*sin(lat_c)^2);

rI = Rs*[cos(lat_c)*cos(long);cos(lat)*sin(long);sin(lat)]; % Target position in ECI
rT = Reci2ned*rI; 

gT = -Mu*rT/(norm(rT))^3;

param.rT = rT;

%% Initial Conditions

x0 = 40*1000;
y0 = 38*1000;
z0 = -120*1000;

vx0 = 0.5;
vy0 = 1;
vz0 = 0.1;

m0 = 47653.153;

ic0 = [x0;y0;z0;vx0;vy0;vz0;m0];
%% FInal conditions

xf = 0;
yf = 0;
zf = 0;
vxf = 0.1;
vyf = 0.5;
vzf = 0.2;

mf = param.mdry;

icf = [xf;yf;zf;vxf;vyf;vzf;mf];

%%


