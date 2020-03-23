clear all; close all; clc

addpath(genpath('MA2_libs'));

cfg = load('TestParam.mat');   % load configFile
params = cfg.TestParam;                    % get the set of parameters

%% Space Angle Duality
load("H_LOS_500.mat");
load("H_NLOS_500.mat");

aLOS = fft(HLOS);
aNLOS = fft(HNLOS);

[n,m] = size(aLOS);

SLOS = aLOS.*conj(aLOS)/(2*pi);
SNLOS = aNLOS.*conj(aNLOS)/(2*pi);

SLOS = mean(SLOS,2);
SNLOS = mean(SNLOS,2);

costheta = [-1 0 1 0];
theta = [-pi -pi/2 0 pi/2];
beta = 2*pi*params.Fc/physconst('LightSpeed');
u = beta*costheta;
u2 = u.*u;

PLOS = sum(SLOS);
PNLOS = sum(SNLOS);

umLOS = sum(u.*SLOS.')./PLOS;
umNLOS = sum(u.*SNLOS.')./PNLOS;

varLOS = sum(u2.*SLOS.')./PLOS;
varNLOS = sum(u2.*SNLOS.')./PNLOS;

sigmaLOS = sqrt(varLOS-(umLOS.*umLOS));
sigmaNLOS = sqrt(varNLOS-(umNLOS.*umNLOS));

thetaLOS=acos(sigmaLOS/beta);
thetaLOS = rad2deg(thetaLOS);

thetaNLOS = acos(sigmaNLOS/beta);
thetaNLOS = rad2deg(thetaNLOS);


RLOS = ifft(SLOS);
RNLOS = ifft(SNLOS);
figure
subplot(2,2,1)
stem(theta,SLOS);
title("LOS Power Angular Spectrum");
xlabel("theta");
xticks([-pi -pi/2 0 pi/2])
xticklabels({'-\pi','-\pi/2','0','\pi/2'})
ylabel("S(u)");
subplot(2,2,2)
stem(theta,SNLOS);
title("NLOS Power Angular Spectrum");
xlabel("theta");
xticks([-pi -pi/2 0 pi/2])
xticklabels({'-\pi','-\pi/2','0','\pi/2'})
ylabel("S(u)");
subplot(2,2,3);
plot(abs(RLOS))
title("LOS Spatial Correlation");
ylabel("|R(\Deltaz)|")
subplot(2,2,4);
plot(abs(RNLOS))
ylabel("|R(\Deltaz)|")
title("NLOS Spatial Correlation");

