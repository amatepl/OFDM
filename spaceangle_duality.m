clear all; close all; clc

addpath(genpath('MA2_libs'));

cfg = load('TestParam.mat');   % load configFile
params = cfg.TestParam;                    % get the set of parameters

%% Space Angle Duality
% load("H_LOS_500.mat");
% load("H_NLOS_500.mat");

load('H_LOS.mat');
HLOS_370 = H(:,370,:);
HLOS_370 = reshape(HLOS_370,size(HLOS_370,1),size(HLOS_370,3));
HLOS_550 = H(:,550,:);
HLOS_550 = reshape(HLOS_550,size(HLOS_550,1),size(HLOS_550,3));
load('H_NLOS.mat')
HNLOS_370 = H(:,370,:);
HNLOS_370 = reshape(HNLOS_370,size(HNLOS_370,1),size(HNLOS_370,3));
HNLOS_550 = H(:,550,:);
HNLOS_550 = reshape(HNLOS_550,size(HNLOS_550,1),size(HNLOS_550,3));

aLOS_370 = fft(HLOS_370);
aLOS_550 = fft(HLOS_550);

aNLOS_370 = fft(HNLOS_370);
aNLOS_550 = fft(HNLOS_550);

[n,m] = size(aLOS_370);

SLOS_370 = aLOS_370.*conj(aLOS_370)/(2*pi);
SLOS_550 = aLOS_550.*conj(aLOS_550)/(2*pi);

SNLOS_370 = aNLOS_370.*conj(aNLOS_370)/(2*pi);
SNLOS_550 = aNLOS_550.*conj(aNLOS_550)/(2*pi);

SLOS_370 = mean(SLOS_370,2);
SLOS_550 = mean(SLOS_550,2);

SNLOS_370 = mean(SNLOS_370,2);
SNLOS_550 = mean(SNLOS_550,2);

costheta = [-1 0 1 0];
% theta=linspace(0,pi,5);
%theta=theta([1 2 3 4]);
theta = [-pi -pi/2 0 pi/2];
beta = 2*pi*params.Fc/physconst('LightSpeed');
u = beta*costheta;
%u = linspace(-0.5,0.5,4);
u2 = u.*u;

PLOS_370 = trapz(SLOS_370);
PLOS_550 = trapz(SLOS_550);

PNLOS_370 = trapz(SNLOS_370);
PNLOS_550 = trapz(SNLOS_550);

umLOS_370 = trapz(u.*SLOS_370.')./PLOS_370;
umLOS_550 = trapz(u.*SLOS_550.')./PLOS_550;

umNLOS_370 = trapz(u.*SNLOS_370.')./PNLOS_370;
umNLOS_550 = trapz(u.*SNLOS_550.')./PNLOS_550;

varLOS_370 = trapz(u2.*SLOS_370.')./PLOS_370;
varLOS_550 = trapz(u2.*SLOS_550.')./PLOS_550;

varNLOS_370 = trapz(u2.*SNLOS_370.')./PNLOS_370;
varNLOS_550 = trapz(u2.*SNLOS_550.')./PNLOS_550;

sigmaLOS_370 = sqrt(varLOS_370-(umLOS_370.*umLOS_370));
sigmaLOS_550 = sqrt(varLOS_550-(umLOS_550.*umLOS_550));

sigmaNLOS_370 = sqrt(varNLOS_370-(umNLOS_370.*umNLOS_370));
sigmaNLOS_550 = sqrt(varNLOS_550-(umNLOS_550.*umNLOS_550));

thetaLOS_370=acos(sigmaLOS_370/beta);
thetaLOS_370 = rad2deg(thetaLOS_370);
thetaLOS_550=acos(sigmaLOS_550/beta);
thetaLOS_550 = rad2deg(thetaLOS_550);

thetaNLOS_370 = acos(sigmaNLOS_370/beta);
thetaNLOS_370 = rad2deg(thetaNLOS_370);
thetaNLOS_550 = acos(sigmaNLOS_550/beta);
thetaNLOS_550 = rad2deg(thetaNLOS_550);


RLOS_370 = ifft(SLOS_370);
RLOS_550 = ifft(SLOS_550);

RNLOS_370 = ifft(SNLOS_370);
RNLOS_550 = ifft(SNLOS_550);

figure
subplot(2,2,1)
stem(theta,SLOS_370);
hold on;
stem(theta,SLOS_550);
legend('subcrr = 370','subcrr = 550');
title("LOS Power Angular Spectrum");
xlabel("theta");
xticks([-pi -pi/2 0 pi/2])
xticklabels({'-\pi','-\pi/2','0','\pi/2'})
ylabel("S(u)");
subplot(2,2,2)
stem(theta,SNLOS_370);
hold on;
stem(theta,SNLOS_550);
legend('subcrr = 370','subcrr = 550');
title("NLOS Power Angular Spectrum");
xlabel("theta");
xticks([-pi -pi/2 0 pi/2])
xticklabels({'-\pi','-\pi/2','0','\pi/2'})
ylabel("S(u)");
subplot(2,2,3);
plot(abs(RLOS_370))
hold on;
plot(abs(RLOS_550));
title("LOS Spatial Correlation");
ylabel("|R(\Deltaz)|")
xlabel("\Deltaz index");
legend('subcrr = 370','subcrr = 550');
subplot(2,2,4);
plot(abs(RNLOS_370))
hold on;
plot(abs(RNLOS_550))
legend('subcrr = 370','subcrr = 550');
ylabel("|R(\Deltaz)|")
xlabel("\Deltaz index");
title("NLOS Spatial Correlation");

