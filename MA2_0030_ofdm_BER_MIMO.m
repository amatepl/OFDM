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

H1 = fftshift(H_NLOS_G1,2);
H1 = H1(:,params.ActiveQIndex);
H1 = permute(H1,[1 3 2]);

H2 = fftshift(H_NLOS_G6,2);
H2 = H2(:,params.ActiveQIndex);
H2 = permute(H2,[1 3 2]);

H = [H1 H2];
W = permute(H,[2 1 3]);

% Compute the pseudo inverse of the H for MIMO
for i = 1:size(W,3)
    W(:,:,i) = pinv(H(:,:,i));
end

dispConfigFile_Test(params);                 % display the parameters

% --- Local parameters ----
Nsymb_ofdm = 10;     % number OFDM symbols to transmit
NsimPerSNR = 50;    % number of simulations per SNR value
Nbits = params.nData * params.nActiveQ * Nbps;
Nr = 4;                             % number of receivers

% define storage variables:
BER_i = zeros(NsimPerSNR,length(SNR_list));

figure, hold on;
plot(abs(H_LOS_G1(:,:).'))
legend('1','2','3','4');
title('Groupe 1');

figure, hold on;
plot(abs(H_NLOS_G1(:,:).'))
legend('1','2','3','4');
title('Groupe 1 - NLOS');


% figure, hold on;
% plot(abs(H_LOS_G2(:,:).'))
% legend('1','2','3','4');
% title('Groupe 2');

figure, hold on;
plot(abs(H_LOS_G6(:,:).'))
legend('1','2','3','4');
title('Groupe 6');

figure, hold on;
plot(abs(H_NLOS_G6(:,:).'))
legend('1','2','3','4');
title('Groupe 6');

progress_indx = 0;
for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = SNR_list(SNR_idx);
            
        % 1. Message, preamble and pilot construction
        [Preamble1, bits_data1, ~] = build_message_test(params,Nbits,Nbps);
        [Preamble2, bits_data2, ~] = build_message_test(params,Nbits,Nbps);

        %bits_tx = vertcat(Preamble,bits_data);
        bits_tx1 = bits_data1;
        bits_tx2 = bits_data2;

        % 2. Modulation of the preamble, message and pilot
        [Qsymb_pre1] = modulation(1,Preamble1,'bpsk');      % Preamble modulation
        [Qsymb_data1] = modulation(Nbps,bits_data1,'qpsk');    % Message modulation
        
        [Qsymb_pre2] = modulation(1,Preamble2,'bpsk');      % Preamble modulation
        [Qsymb_data2] = modulation(Nbps,bits_data2,'qpsk');    % Message modulation
        
        % OFDM symbols [2 x preamble + message,1]
        Qsymb_tx1 = vertcat(Qsymb_pre1,Qsymb_data1);  
        Qsymb_tx2 = vertcat(Qsymb_pre2,Qsymb_data2);

        % 3. OFDM Transmitter: 
        [signal_tx1] = transmitter4(params, Qsymb_pre1, Qsymb_data1);
        preamble1 = Qsymb_pre1(1:params.nActiveQ);    % Take only one copy of the preamble
        
        [signal_tx2] = transmitter4(params, Qsymb_pre2, Qsymb_data2);
        preamble2 = Qsymb_pre2(1:params.nActiveQ);    % Take only one copy of the preamble
        
        % 4. Channel propagation: 
          
%         signal_rx_1 = channel_propagation4(params,signal_tx1,H_LOS_G1(1:4,:),SNR,Nr);
%         signal_rx_2 = channel_propagation4(params,signal_tx2,H_LOS_G6(1:Nr,:),SNR,Nr);
        
        signal_rx = channel_propagationMIMO(params,signal_tx1,signal_tx2,H_NLOS_G1(1:4,:),H_NLOS_G6(1:Nr,:),SNR,Nr);
        
%         figure, hold on;
%         plot(abs(signal_rx_1(:,1,1)));
%         title('signal rx1');
%         
%         figure, hold on;
%         plot(abs(signal_rx_2(:,1,1)));
%         title('signal rx2');
        
%         signal_rx = signal_rx_1+signal_rx_2;
%         
%         figure, hold on;
%         plot(abs(signal_rx(:,1,1)));
        
        [Qsymb_rx_1,Qsymb_rx_2] = receiverMIMO(params,signal_rx,W,H_NLOS_G1(1:4,:),H_NLOS_G6(1:Nr,:) ,Nr);
        
        % USER 1
        % 5. Demodulation:
        bits_rx_1 = demodulation(params,Qsymb_rx_1(2*params.nActiveQ+1:end),'qpsk');
        % compute BER
        bitErrorRate_1 = sum(abs(bits_tx1 - bits_rx_1),'all');
        
        % USER 2
       
        % 5. Demodulation:
        bits_rx_2 = demodulation(params,Qsymb_rx_2(2*params.nActiveQ+1:end),'qpsk');       
        % compute BER
        bitErrorRate_2 = sum(abs(bits_tx2 - bits_rx_2),'all');
        
        
          
%% Initial code        
%         % 1. QAM Modulation.
%         [bits_tx,Qsymb_tx] = modulation(params,Nbits);
% 
%         % 2. OFDM Transmitter: 
%         signal_tx = transmitter(params,Qsymb_tx,Nsymb_ofdm);
% 
%         % 3. Channel propagation: 
%         signal_rx = channel_propagation(params,signal_tx,SNR);
% 
%         % 4. OFDM Receiver:
%         Qsymb_rx = receiver(params,signal_rx,Nsymb_ofdm);
% 
%         % 5. Demodulation:
%         bits_rx = demodulation(params,Qsymb_rx);
%%
        
        
        BER_i(sim_idx,SNR_idx) = (bitErrorRate_1+bitErrorRate_2 )/ (2*length(bits_tx1));
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
xlim([-5 15]);