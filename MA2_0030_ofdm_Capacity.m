% Opera - ULB
% February 2020
%
% This script generates the BER curve using the same functions used as in 
% the script: 'MA2_0020_ofdm_structure.m', Please be sure that your code is
% running well in the previous script before starting this part.
% -------------------------------------------------------------------------
clear; close all; clc;
addpath(genpath('MA2_libs'));           % add libraries
addpath(genpath('../Channel_frequency_response_LOSNLOS-each_group')); % add chanel responses of different groups

cfg = load('MA2_lab_parameters.mat');   % load configFile
SNR_list = cfg.params.SNR_list;                    % get the set of parameters
Nbps = cfg.params.modulation.Nbps; 
cfg = load('TestParam4.mat');   % load configFile
params = cfg.TestParam4;                    % get the set of parameters
H_LOS_G1 = load('H_LOS_G1.mat');
H_LOS_G1 = H_LOS_G1.H;
H_NLOS_G1 = load('H_NLOS_G1.mat');
H_NLOS_G1 = H_NLOS_G1.H;
H_LOS_G2 = load('H_LOS_G2.mat');
H_LOS_G2 = H_LOS_G2.H_los_g2;
H_LOS_G6 = load('Hest_LOS.mat');
H_LOS_G6 = H_LOS_G6.Hest;
H_LOS_G6 = H_LOS_G6.';
z = zeros(4,params.Q);
z(:,params.ActiveQIndex) = (H_LOS_G6);
H_LOS_G6 = ifftshift(z);
H_NLOS_G6 = load('Hest_NLOS.mat');
H_NLOS_G6 = H_NLOS_G6.Hest;
H_NLOS_G6 = H_NLOS_G6.';
z = zeros(4,params.Q);
z(:,params.ActiveQIndex) = (H_NLOS_G6);
H_NLOS_G6 = ifftshift(z);

clear z;

Htype = 'NLOS';                     % NLOS or LOS

switch (Htype)
    case 'NLOS'
        H1 = H_NLOS_G1;
        H2 = H_NLOS_G6;
    case 'LOS'
        H1 = H_LOS_G1;
        H2 = H_LOS_G6;
end

H1 = fftshift(H1,2);
H1 = H1(:,params.ActiveQIndex);
H1 = permute(H1,[1 3 2]);

H2 = fftshift(H2,2);
H2 = H2(:,params.ActiveQIndex);
H2 = permute(H2,[1 3 2]);

H = [H1 H2];
W = permute(H,[2 1 3]);
[U,S,V] = svd(H1(1:4,:,5));
% Compute the pseudo inverse of the H for MIMO
for i = 1:size(W,3)
    W(:,:,i) = pinv(H(:,:,i));
end

dispConfigFile_Test(params);                 % display the parameters

% --- Local parameters ----

Nbits = params.nData * params.nActiveQ * Nbps;
Nr = 4;                             % number of receivers
Nt = 1;                             % number of transmitters

% SIMO

C = zeros(4,length(SNR_list));

progress_indx = 0;
for nAntenas = 1:4
    progress_indx =  progress_indx + 1;
    alpha = zeros(1,params.nActiveQ);
    for Qi = 1:params.nActiveQ
        [~,S,~] = svd(H1(1:nAntenas,:,Qi));
        alpha(Qi) = S(1);
    end
    for SNR_idx = 1:length(SNR_list)
        SNR = SNR_list(SNR_idx);
        C(nAntenas,SNR_idx) = mean(log2(1 + (10^(SNR/10)*alpha.^2)));
    end
    
    progress = round(progress_indx/4*100);
    disp([' Simulation executed @ ', num2str(progress), '%']);
    disp('');
end

% MIMO

C_MIMO = zeros(1,length(SNR_list));
alpha = zeros(2,params.nActiveQ);
    for Qi = 1:params.nActiveQ
        [~,S,~] = svd(H(1:4,:,Qi));
        alpha(:,Qi) = S(1:2);
    end
    for SNR_idx = 1:length(SNR_list)
        SNR = SNR_list(SNR_idx);
        C_MIMO(1,SNR_idx) = mean(sum(log2(1 + (10^(SNR/10)*alpha.^2)/2),1));
    end

% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------

disp('$$ Displaying results:');
figure('Name',join(['Capacity',Htype]));hold on;
% semilogy(SNR_list,mean(BER_i,1));
plot(SNR_list,C.');
plot(SNR_list,C_MIMO.')
grid on; hold on;
% ber theoretical
ber_theo = berawgn(SNR_list,'qam',2^(Nbps));
% semilogy(SNR_list,ber_theo,'--');
legend('M_T = 1, M_R = 1','M_T = 1, M_R = 2','M_T = 1, M_R = 3','M_T = 1, M_R = 4','M_T = 2, M_R = 4');
xlabel('SNR dB');ylabel('Capacity in bps/Hz');
xlim([-5 15]);
title(join(['Shanon capacity for ',Htype,' propagation']));