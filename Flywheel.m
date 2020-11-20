%% Mini-project: Conversion chain design for a flywheel

%% Default commands
close all; clear all; clc;

%% Load data
load("Centrale250kW_16_sept.mat");
Ppv = Grandeur_lapalud_16sept(:, 4);
Fe = 0.0069; % sampling frequency
Te = 1/Fe;
n = 1440;
t = 0 : Te : (n - 1)*Te;
coeff = 0.95;

%M = round((3.3/12)-1)
m = 477;
B = fir1(m, 2*Fe, "low", hamming(m+1)); % design a low-pass filter
PV_filtre = filtfilt(B, 1, Ppv);

k = 0.5;
Sigmax = 1800*10^6;
ro = 7800;

Pmoy_nf = mean(Ppv); %Puissance moyenne est la même entre les
Pmoy_f = mean(PV_filtre);
Pmax_nf = max(Ppv);
Pmax_f = max(PV_filtre);

PHP = ones(n,1);
Pmax_PHP = zeros(n,1);
lb = 0;
ub = 1;

for i = 600 : 900
   Pmax_PHP(i, 1) = abs(Ppv(i, 1) - PV_filtre(i, 1));
   if PV_filtre(i, 1) >= 0
       PHP(i, 1) = 1 - PV_filtre(i, 1)/Pmax_PHP(i, 1);
   else
       PHP(i, 1) = 1;
   end  
end

P_supplied = zeros(n, 1); 
P_stored = zeros(n, 1);

for i = 1 : n
    if Ppv(i, 1) - PV_filtre(i, 1) < 0
        P_supplied(i, 1) = (Ppv(i, 1) - PV_filtre(i, 1))/coeff;
    else
        P_stored(i, 1) = (Ppv(i, 1) - PV_filtre(i, 1))*coeff;
    end
end

P = P_supplied + P_stored;
E = zeros(n, 1);
E0 = 0;
for i = 1 : n
    E(i, 1) = E0 - trapz(P(1: i));
end

Eu = (max(E) - min(E));
Eu_units = (max(E) - min(E))*60;

E0p = (0.3/0.7)*Eu - min(E);
Es = zeros(n, 1);
for i = 1 : n
    Es(i, 1) = E0p - trapz(P(1: i));
end

Min_Cap = zeros(n, 1);
for i = 1 : n
    Min_Cap(i, 1) = (0.3/0.7)*Eu;
end

Es_unit = Es*60;

Emin = 0; % corresponding to a zero speed

DoDmax = 0.7;
Soc = 1 - DoDmax;
capacity = Eu_units/DoDmax;
w_min = 2760*2*pi/60;
w_max = w_min/sqrt(Soc);
J=2*capacity/(w_max^2);

%%
L = 0.15; % length of the motor
R = (J/(ro*k*pi*L))^(1/4); % radius
M = ro*pi*R*R*L; % weight

%%
%Graph that represents the production
figure;
grid on;
hold on;
plot(t, Ppv);
plot(t, PV_filtre);
plot(t, P_supplied);
plot(t, P_stored);
title('Different types of production');
xlabel('Time (s)');
ylabel('Production (W)');
legend('Production not filtered', 'Production filtered', 'Production supplied', 'Production stored');
hold off;

% PHP
figure;
grid on;
hold on;
title('PHP and Pmax');
plot(t, PV_filtre);
plot(t, abs(PHP));
xlabel('Time (s)');
ylabel('PHP');
legend('PHP','Pmax');
hold off;

% Graph after the capacity
figure;
grid on;
hold on;
title('Shifting the energy stored');
plot(E);
plot(Es);
plot(Min_Cap)
xlabel('Time (s)');
ylabel('Energy');
legend('Energy calculated', 'Energy after shift (higher than 30% capacity)', 'Line of 30% capacity');
hold off;

vect_speed = sqrt(2*Es_unit/J);
vect_speed_unit = vect_speed *(60/(2*pi)); % round/min

Torque_f = P./vect_speed;

Torque_max1 = zeros(n, 1);
for i = 1 : n
    Torque_max1(i) = 240; % maximum torque
end
Torque_max4 = -Torque_max1;

Torque_max2 = zeros(n, 1);
for i = 1 : n
    Torque_max2(i) = max(P)/vect_speed(i);
end

Torque_max3 = -Torque_max2;

figure;
hold on;
scatter(vect_speed_unit, Torque_f);
plot(vect_speed_unit, Torque_max1);
plot(vect_speed_unit, Torque_max2);
plot(vect_speed_unit, Torque_max3);
plot(vect_speed_unit, Torque_max4);
grid on;
legend('Torque/speed', 'Limit line 1', 'Limit line 2', 'Limit line 3', 'Limit line 4');
xlabel('Speed (rpm)');
ylabel('Torque (N.m)');
title('Torque in function of speed');
hold off;

figure;
plot(t, vect_speed_unit);
grid on;
xlabel('Time (s)');
ylabel('Speed (rad/s)');
title('Speed in function of time');
legend('Speed of the motor');

Torque_supplied = P_supplied./vect_speed;
Torque_stored = P_stored./vect_speed;
Torque_moteur = Torque_supplied + Torque_stored;

% Motor torque
figure;
plot(t, Torque_moteur);
grid on;
xlabel('Time (S)');
ylabel('Torque (N.m)');
title('Torque in function of time');
legend('Torque of the motor');

% Calculate isq, isd. We take delta = pi/2
delta = pi/2;
p = 4;
ku = 1.17;
Fifsd = ku/p;
isq = Torque_moteur./(p*Fifsd);
isd = 0*isq;
iseff = sqrt((isq.*isq + isd.*isd)/3);

%Calculate the voltages vsd, vsq
Ld = 0.00259;
Lq = 0.00262;
Rs = 0.0952;
ws = p*vect_speed;
vsd = Rs*isd - ws*Lq.*isq;
vsq = Rs*isq + ws*Ld.*isd + ws*Fifsd;
vseff = sqrt((vsq.*vsq + vsd.*vsd)/3);

%% Plot isd,isq and vsd,vsq :
figure;
subplot(211);
hold on;
plot(t, isq);
plot(t, isd);
grid on;
legend('isq','isd');
xlabel('Time (s)');
ylabel('Current (A)');
title('Current in function of time');
hold off;
subplot(212);
hold on;
plot(t, vsd);
plot(t, vsq);
grid on;
legend('vsd','vsq');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('Voltage in function of time');
hold off;

%% Puissance absorbée
Pabs = vsd.*isd + vsq.*isq;

figure;
plot(t, Pabs);
grid on;
xlabel('Time (s)');
ylabel('Power stored (W)');
title('Power stored in function of time');
legend('Power stored');

%% Efficiency
efficiency = zeros(n, 1);
for i = 1 : n
    if Pabs(i) < 0
        efficiency(i) = abs(Pabs(i))/abs(P(i));
    else
        efficiency(i) = abs(P(i))/abs(Pabs(i));
    end
end

figure;
plot(t, efficiency);
grid on;
xlabel('Time (s)');
ylabel('Efficiency');
title('Efficiency in function of time');
legend('Efficiency');

% Control : we are trying to control the currents isd and isq
% Indeed, by doing this, we control vsd and vsq at the same time
% Open-loop transfer function (TF):
Gd = tf([1/Rs], [Ld/Rs 1]); % TF for isd
Gq = tf([1/Rs], [Lq/Rs 1]); % TF for isq

% PI controller for the current
ksides = 0.7; 

trdes_d = (1/10)*3*Ld/Rs; % We find that with this value, this settling time (st) is 4 times greater than st min and 10x smaller than st max (the st in open loop)
wn_d = 4/trdes_d;
Kpd = 2*ksides*Ld*wn_d - Rs;
Tid = Kpd/(Ld*wn_d*wn_d);
PId = tf([Kpd*Tid Kpd], [Tid 0]);

trdes_q = (1/10)*3*Ld/Rs;
wn_q = 4/trdes_q;
Kpq = 2*ksides*Lq*wn_q - Rs;
Tiq = Kpq/(Lq*wn_q*wn_q);
PIq = tf([Kpq*Tiq Kpq], [Tiq 0]);

% Closed-loop system
Hbfd = zpk(PId*Gd/(1 + PId*Gd));
Hbfq = zpk(PIq*Gq/(1 + PIq*Gq));

figure;
grid on;
subplot(221);
step(Gd);
legend('Gd : transfer function, open loop');
subplot(222);
step(Gq);
legend('Gq : transfer function, open loop');
subplot(223);
step(Hbfd);
legend('Hbfd : transfer function, close loop');
subplot(224);
step(Hbfq);
legend('Hbfq : transfer function, close loop');

% Speed control
vect_speed_d = zeros(n, 1);
Te1 = 60;
for i = 1 : (n - 1)
    vect_speed_d(i) = (vect_speed(i + 1) - vect_speed(i))/Te1;
end

omega_bis = [ones(n, 1) vect_speed]; % We find torque in the form: torque = f0 + f1*vec_speed -> create a matrix to use least squares identification
f = omega_bis\(Torque_f - J*vect_speed_d);
f0 = f(1);
f1 = f(2);

figure;
grid on;
hold on;
scatter(vect_speed, Torque_f);
plot(vect_speed, omega_bis*f);
title('Approximation of the relation between the torque and the speed');
legend('Real graph', 'Least squares approximation');
xlabel('Speed');
ylabel('Torque');
hold off;

G_omega = tf([1/f1], [J/f1 1]); % TF of system for omega

%% PI controller for the speed

trdes_omega = 10;
wn_omega = 4/trdes_omega;
Kp_omega = 2*ksides*J*wn_omega - f1;
Ti_omega = Kp_omega/(J*wn_omega*wn_omega);
PI_omega = tf([Kp_omega*Ti_omega Kp_omega], [Ti_omega 0]);

% tau_omega = trdes_omega/3;
% Kp_omega = 1;
% Ti_omega = tau_omega/2;
% PI_omega = tf( [Kp_omega*Ti_omega Kp_omega], [Ti_omega 0]);

% Close-loop system
Hbf_omega = zpk(PI_omega*G_omega/(1 + PI_omega*G_omega));

figure;
grid on;
subplot(211);
step(G_omega);
legend('Gd : transfer function, open loop');
subplot(212);
step(Hbf_omega);
legend('Hbf_omega : transfer function, close loop');

%%
m = 60;
k = 100;
reference_speed = zeros(n + m + k, 2);
for i = 1 : size(reference_speed)
    reference_speed(i, 1) = i*Te;
end
for i = 1 : m
    reference_speed(i, 2) = 6*i;
end
for i = (m + 1) : (m + n)
    reference_speed(i, 2) = vect_speed(i - m);
end

%%
k = 60;
reference_speed1 = zeros(n + k, 2);
for i = 1 : k
        reference_speed1(i, 1) = (i);
        reference_speed1(i, 2) = 6*i;
end

for i = k + 1 : size(reference_speed1)
        reference_speed1(i, 1) = (i);
        reference_speed1(i, 2) = vect_speed(i - k);
end

%% Attempt to stop
m = 100;
reference_speed2 = zeros(n + m, 2);
for i = 1 : (m + n)
    reference_speed2(i, 1) = (i);
end
for i = 1 : n
    reference_speed2(i,2) = vect_speed(i);
end

%% Comparison
densityPower = P_stored/M;
figure;
plot(densityPower);
grid on;
legend('Density power');
title('Density power in function of time');
xlabel('Time (s)');
ylabel('Density power (W/kg)');