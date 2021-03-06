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

dispConfigFile_Test(params);                 % display the parameters

% --- Local parameters ----

%% Parameters to set ---------------------------------------
NsimPerSNR = 10;    % number of simulations per SNR value
Nr = 4;                             % number of receivers
Htype = 'NLOS';                      % NLOS or LOS
%% ---------------------------------------------------------
% params.nData = 50;

Nsymb_ofdm = 10;     % number OFDM symbols to transmit
Nbits = params.nData * params.nActiveQ * Nbps;

params.Nbps = 2;    % Nbps
params.modulation = 'qpsk';

switch (Htype)
    case 'NLOS'
        H = H_NLOS_G1;
    case 'LOS'
        H = H_LOS_G1;
end


Hw = fftshift(H,2);
Hw = Hw(:,params.ActiveQIndex);
Hw = permute(Hw,[1 3 2]);

W = permute(Hw,[2 1 3]);

% Compute the pseudo inverse of the H for MIMO
for i = 1:size(W,3)
    W(:,:,i) = pinv(Hw(:,:,i));
end

% define storage variables:
BER_i = zeros(NsimPerSNR,length(SNR_list));

figure, hold on;
plot(abs(H(:,:).'))
legend('1','2','3','4');
title('Groupe 1');

% figure, hold on;
% plot(abs(H_NLOS_G1(:,:).'))
% legend('1','2','3','4');
% title('Groupe 1 - NLOS');



progress_indx = 0;
for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = SNR_list(SNR_idx);
            
        % 1. Message, preamble and pilot construction
        [Preamble, bits_data, ~] = build_message_test(params,Nbits);

        %bits_tx = vertcat(Preamble,bits_data);
        bits_tx = bits_data;

        % 2. Modulation of the preamble, message and pilot
        [Qsymb_pre] = modulation(params,Preamble,'qpsk');      % Preamble modulation
        [Qsymb_data] = modulation(params,bits_data,'qpsk');    % Message modulation
        
        % OFDM symbols [2 x preamble + message,1]
        Qsymb_tx = vertcat(Qsymb_pre,Qsymb_data);             

        % 3. OFDM Transmitter: 
        [signal_tx] = transmitter4(params, Qsymb_pre, Qsymb_data);
        preamble = Qsymb_pre(1:params.nActiveQ);    % Take only one copy of the preamble
        
        % 4. Channel propagation: 
        
        % USER 1
        signal_rx_1 = channel_propagation4(params,signal_tx,H(1:Nr,:),SNR,Nr);
        [hz,Qsymb_rx_1] = receiver4(params,signal_rx_1,params.nData, preamble,Nr,W);
        % 5. Demodulation:
        bits_rx_1 = demodulation(params,Qsymb_rx_1(2*params.nActiveQ+1:end),'qpsk');
        % compute BER
        bitErrorRate_1 = sum(abs(bits_tx - bits_rx_1),'all');
        
        BER_i(sim_idx,SNR_idx) = (bitErrorRate_1)/(length(bits_tx));
    end
    progress = round(progress_indx/NsimPerSNR/length(SNR_list)*100);
    disp([' Simulation executed @ ', num2str(progress), '%']);
    disp('');
end

% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------

disp('$$ Displaying results:');
figure;
semilogy(SNR_list,mean(BER_i,1));
grid on; hold on;
% ber theoretical
ber_theo = berawgn(SNR_list,'qam',2^(Nbps));
semilogy(SNR_list,ber_theo,'--');
legend('myBER','theoretical');
xlabel('SNR dB');ylabel('Probability of error');
xlim([-5 20]);
ylim([10^(-6) 1]);

fpath = './Results';
filename = join(['SIMO_',Htype,'_',num2str(Nr),'.mat']);
save(fullfile(fpath,Htype,filename), 'BER_i','-mat');