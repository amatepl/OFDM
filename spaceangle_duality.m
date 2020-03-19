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

costheta = [1 0 -1 0];
beta = 2*pi*params.Fc/physconst('LightSpeed');
u = beta*costheta;

PLOS = sum(SLOS);
PNLOS = sum(SNLOS);

u = repmat(u.',[1,m]);
u2 = u.*u;

umLOS = sum(u.*SLOS)./PLOS;
umNLOS = sum(u.*SNLOS)./PNLOS;

varLOS = sum(u2.*SLOS)./PLOS;
varNLOS = sum(u2.*SNLOS)./PNLOS;

sigmaLOS = sqrt(varLOS-(umLOS.*umLOS));
sigmaNLOS = sqrt(varNLOS-(umNLOS.*umNLOS));

thetaLOS=acos(sigmaLOS/beta);
thetaLOS = rad2deg(thetaLOS);

thetaNLOS = acos(sigmaNLOS/beta);
thetaNLOS = rad2deg(thetaNLOS);

RLOS = ifft(SLOS);
RNLOS = ifft(SNLOS);
figure
subplot(1,2,1);
plot(abs(RLOS(:,1:10)))
title("LOS Spatial Correlation");
subplot(1,2,2);
plot(abs(RNLOS(:,1:10)))
title("NLOS Spatial Correlation");

