clear all
clc

%% VCO

Kvco = 5e+6; % Hz/V
fo = 1000; % Hz

%% Loop filter

Kp = 800;
Ki = 80;
fc = 90;
Wc = 2*pi*fc;
tau = 1/Wc;
%% 
t1 = tf([Kp Ki],[1 Wc 0])
bp = bodeplot(t1)
bp.FrequencyUnit = "Hz";

%% Digital PFD
% Icp = 50e-6; %uA
% Tup = 100e-9; %ns
% C = 1e-6; %uF
Kpd_scale = 1;

%% Input
fref = 500; %Hz