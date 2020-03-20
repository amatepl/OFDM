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
H_LOS_G2 = load('H_LOS_G2.mat');
H_LOS_G2 = H_LOS_G2.H_los_g2;
H_LOS_G6 = load('Hest_LOS.mat');
H_LOS_G6 = H_LOS_G6.Hest;
H_LOS_G6 = H_LOS_G6.';
z = zeros(4,params.Q);
z(:,params.ActiveQIndex) = H_LOS_G6;
H_LOS_G6 = z;
clear z;

dispConfigFile_Test(params);                 % display the parameters

% --- Local parameters ----
Nsymb_ofdm = 10;     % number OFDM symbols to transmit
NsimPerSNR = 10;    % number of simulations per SNR value
Nbits = params.nData * params.nActiveQ * Nbps;
Nr = 4;                             % number of receivers

% define storage variables:
BER_i = zeros(NsimPerSNR,length(SNR_list));

figure, hold on;
plot(abs(H_LOS_G1(1,:)))

progress_indx = 0;
for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = SNR_list(SNR_idx);
            
        % 1. Message, preamble and pilot construction
        [Preamble, bits_data, ~] = build_message_test(params,Nbits,Nbps);
%         Preamble = ones(params.nActiveQ*2,1);

        %bits_tx = vertcat(Preamble,bits_data);
        bits_tx = bits_data;

        % 2. Modulation of the preamble, message and pilot
        [Qsymb_pre] = modulation(1,Preamble,'bpsk');      % Preamble modulation
        [Qsymb_data] = modulation(Nbps,bits_data,'qpsk');    % Message modulation
        
        % OFDM symbols [2 x preamble + message,1]
        Qsymb_tx = vertcat(Qsymb_pre,Qsymb_data);             

        % 3. OFDM Transmitter: 
        [signal_tx] = transmitter4(params, Qsymb_pre, Qsymb_data);
%         preambleLCP = signal_tx(:,1:params.Q+params.LCP);

        % 4. Channel propagation: 
        signal_rx = channel_propagation4(params,signal_tx,H_LOS_G1(1:Nr,:),SNR,Nr);
%         y = testFunc(signal_tx,H_LOS_G1);
        
        %Average over the antennas
        %STO_estimated = round(mean(STO_estimated,'all'));
        %CFO_estimated = mean(CFO_estimated,'all');

        preamble = Qsymb_pre(1:params.nActiveQ);    % Take only one copy of the preamble

        [hz,Qsymb_rx] = receiver4(params,signal_rx,params.nData, preamble);

        % 5. Demodulation:
        bits_rx = demodulation(params,Qsymb_rx(2*params.nActiveQ+1:end),'qpsk');
        
          
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
        
        % compute BER
        bitErrorRate = sum(abs(bits_tx - bits_rx),'all');
        BER_i(sim_idx,SNR_idx) = bitErrorRate / length(bits_tx);
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