% Opera - ULB
% February 2020
%
% This script generates the BER curve using the same functions used as in 
% the script: 'MA2_0020_ofdm_structure.m', Please be sure that your code is
% running well in the previous script before starting this part.
% -------------------------------------------------------------------------
clear; close all; clc;
addpath(genpath('MA2_libs'));           % add libraries
cfg = load('MA2_lab_parameters.mat');   % load configFile
params = cfg.params;                    % get the set of parameters
dispConfigFile(params);                 % display the parameters

% --- Local parameters ----
Nsymb_ofdm = 10;     % number OFDM symbols to transmit
NsimPerSNR = 10;    % number of simulations per SNR value
Nbits = Nsymb_ofdm * params.ofdm.N_subcrr * params.modulation.Nbps;


% define storage variables:
BER_i = zeros(NsimPerSNR,length(params.SNR_list));


progress_indx = 0;
for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(params.SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = params.SNR_list(SNR_idx);
        % 1. QAM Modulation.
        [bits_tx,Qsymb_tx] = modulation(params,Nbits);

        % 2. OFDM Transmitter: 
        signal_tx = transmitter(params,Qsymb_tx,Nsymb_ofdm);

        % 3. Channel propagation: 
        signal_rx = channel_propagation(params,signal_tx,SNR);

        % 4. OFDM Receiver:
        Qsymb_rx = receiver(params,signal_rx,Nsymb_ofdm);

        % 5. Demodulation:
        bits_rx = demodulation(params,Qsymb_rx);
        
        % compute BER
        bitErrorRate = sum(abs(bits_tx - bits_rx),'all');
        BER_i(sim_idx,SNR_idx) = bitErrorRate / length(bits_tx);
    end
    progress = round(progress_indx/NsimPerSNR/length(params.SNR_list)*100);
    disp([' Simulation executed @ ', num2str(progress), '%']);
    disp('');
end

% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------

disp('$$ Displaying results:');
figure;
semilogy(params.SNR_list,mean(BER_i,1));
grid on; hold on;
% ber theoretical
ber_theo = berawgn(params.SNR_list,'qam',2^(params.modulation.Nbps));
semilogy(params.SNR_list,ber_theo,'--');
legend('myBER','theoretical');
xlabel('SNR dB');ylabel('Probability of error');
xlim([-5 15]);