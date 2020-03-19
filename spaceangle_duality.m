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

theta = [0 pi/2 pi 3*pi/2];
beta = 2*pi*params.Fc/physconst('LightSpeed');
u = beta*cos(theta);

PLOS = sum(SLOS,1);
PNLOS = sum(SNLOS,1);

u = repmat(u.',[1,m]);
u2 = u.*u;
umLOS = sum(u.*SLOS,1)./PLOS;
umNLOS = sum(u.*SNLOS,1)./PNLOS;

varLOS = sum(u2.*SLOS,1)./PLOS;
varNLOS = sum(u2.*SNLOS,1)./PNLOS;

sigmaLOS = sqrt(varLOS-(umLOS.*umLOS));
sigmaNLOS = sqrt(varNLOS-(umNLOS.*umNLOS));

RLOS = ifft(SLOS)/(2*pi);
RNLOS = ifft(SNLOS)/(2*pi);
figure
subplot(1,2,1);
plot(abs(RLOS(:,1:10)./abs(RLOS(1,1:10))))
title("LOS Spatial Correlation");
subplot(1,2,2);
plot(abs(RNLOS(:,1:10)./abs(RNLOS(1,1:10))))
title("NLOS Spatial Correlation");

