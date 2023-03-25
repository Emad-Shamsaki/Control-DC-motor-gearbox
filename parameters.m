%% General parameters and conversion gains
 clear all
 close all
 clc
 % conversion gains
 rpm2rads = 2*pi/60; % [rpm] −> [rad/s]
 rads2rpm = 60/2/pi; % [rad/s] −> [rpm]
 ozin2Nm = 0.706e-2; % [oz*inch] −> [N*m]
 rad2deg = 180/pi; % [rad] −> [deg]
 deg2rad = pi/180;
 degs2rpm = 60/360; % [deg/s] −> [rpm]
%% DC motor nominal parameters
 mot.R = 2.6; % armature resistance
 mot.L = 180e-6; % armature inductance
 mot.Kt = 1.088 * ozin2Nm; % torque constant
 mot.Ke = 0.804e-3 * rads2rpm; % back−EMF constant
 mot.J = 5.523e-5 * ozin2Nm; % rotor inertia
 mot.B = 0.0; % viscous friction coeff(n.a.)
 mot.eta = 0.69; % motor efficiency
 mot.PN = 3.23/mot.eta; % nominal output power
 mot.UN = 6; % nominal voltage
 mot.IN = mot.PN/mot.UN; % nominal current
 mot.tauN = mot.Kt*mot.IN; % nominal torque
 mot.taus = 2.42 * ozin2Nm; % stall torque
 mot.w0 = 7200 * rpm2rads; % no−load speed

  %% Gearbox nominal parameters
 % planetary gearbox Micromotor SA 23/1
 gbox.N = 14; % 1st reduction ratio (planetary gearbox)
 gbox.eta = 0.80; % gearbox efficiency

 %% Mechanical load nominal parameters
 % inertia disc params
 mld.Jl = 3e-5; % load disc inertia
 mld.Bl = 0; % load viscous coeff (n.a.)
 mld.tausf = 7e-3;

 %% shunt resistor
 Rs = 0.5;


%% Voltage driver nominal parameters

 % op−amp circuit params
 drv.R1 = 7.5e3; % op−amp input resistor (dac to non−inverting in)
 drv.R2 = 1.6e3; % op−amp input resistor (non−inverting in to gnd)
 drv.R3 = 1.2e3; % op−amp feedback resistor (output to inverting in)
 drv.R4 = 0.5e3; % op−amp feedback resistor (inverting in to gnd)
 drv.C1 = 100e-9; % op−amp input capacitor
 drv.outmax = 12; % op−amp max output voltage

 % voltage driver dc−gain
 drv.dcgain = drv.R2/(drv.R1+drv.R2) * (1 + drv.R3/drv.R4);

 % voltage driver time constant
 drv.Tc = drv.C1 * drv.R1*drv.R2/(drv.R1+drv.R2);
%%
Req = mot.R + Rs;
Jeq = mot.J + mld.Jl / (gbox.N ^ 2);
Beq = 2.0e-6;
%% PID parameters
Mp = 0.1;
Ts_5 = 0.15;

B_PID = 0;
damp_factor = log(1/Mp)/sqrt(pi^2 + log(1/Mp) ^ 2);
w_gc = 3 / (damp_factor * Ts_5);
phi_m = atan(2 * damp_factor / sqrt(sqrt(1 + 4 * damp_factor ^4) - 2 * damp_factor ^ 2));

Km = drv.dcgain * mot.Kt / (Req * B_PID + mot.Kt * mot.Ke);
Tm = Req * Jeq / (Req * B_PID + mot.Kt * mot.Ke);
plant = tf(Km,gbox.N*[Tm 1 0]);

[mag,phase,w_gc] = bode(plant,w_gc);
deltaK = 1/mag;
delta_phi = - pi + phi_m - deg2rad * phase;

Kp = deltaK * cos(delta_phi) * 2;
Td = (tan(delta_phi)+sqrt((tan(delta_phi))^2 + 1))/ 2 / w_gc;
Ti = 4 * Td;
Kd = Kp * Td / 1.5;
Ki = Kp / Ti / 30;
Tl = 1 /(3.5 * w_gc);

%%
% robust controt
%[A,B,C,D] = tf2ss(Km,gbox.N*[Tm 1 0]);
A = [0 1;0 -1/Tm];
B = [0;Km/(gbox.N*Tm)];
C = [1 0];
D = 0;
A_new = [0,C;[0;0],A];
B_new = [0;B];
C_new = [0,C];

desired_poles = [-2*damp_factor*w_gc+j*w_gc*sqrt(1-damp_factor^2),-2*damp_factor*w_gc-j*w_gc*sqrt(1-damp_factor^2),-2*damp_factor*w_gc];
%P2 = [-zita*w_gc+w_gc*sqrt(1-zita^2)*j -zita*w_gc-w_gc*sqrt(1-zita^2)*j -zita*w_gc];
%P2 = [-zita*w_gc -zita*w_gc -zita*w_gc];
%P2 = [-2*zita*w_gc+w_gc*sqrt(1-zita^2)*j -2*zita*w_gc-w_gc*sqrt(1-zita^2)*j -2*zita*w_gc];
%P2 = [-2*zita*w_gc+w_gc*sqrt(1-zita^2)*j -2*zita*w_gc-w_gc*sqrt(1-zita^2)*j -3*zita*w_gc];
K = place(A_new,B_new,desired_poles);

K1 = [K(1,2),K(1,3)];
K2 = K(1,1);
%%
%Derivative block is a butterworth filter
w_c = 2*pi*50;
sci = 1/sqrt(2);
