% Opera - ULB
% February 2020
%
% This script generates the BER curve using the same functions used as in 
% the script: 'MA2_0020_ofdm_structure.m', Please be sure that your code is
% running well in the previous script before starting this part.
% -------------------------------------------------------------------------
clear ; close all; clc;
addpath(genpath('MA2_libs'));           % add libraries

cfg = load('TestParam.mat');            % load configFile
params = cfg.TestParam;                 % get the set of parameters
dispConfigFile_Test(params);            % display the parameters
params.N_pilots = 126;                  % add the number of pilots
params.N_zeros = 0;                     % add the number of zero ofdm symbol
params.SNR_list = -5:5:15;              % SNR list for BER curve
params.Nbps = 1;                        % Modulation order
params.modulation = 'bpsk';
NsimPerSNR = 100;                        % number of simulations per SNR value

%% --- Local parameters
STO = 0;                                % Time offset (switching unit vector)
% delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
CFO = 0;                                % Carrier frequency offset in ppm
Nr = 1;                                 % number of receivers

% Number of bits knowing the inactive subcarriers and the number of pilots
Nbits = params.nData * (params.nActiveQ-params.N_pilots) * params.Nbps;
frame_size = (params.nPreamble+params.nData+params.N_zeros)*(params.Q+params.LCP);
Nsymb = (params.Q+params.LCP)*(params.nData+params.nPreamble);


%% ------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------

% 1. Message, preamble and pilot construction
[Preamble, bits_data, bits_pilot] = build_message_test(params,Nbits);
 
% 2. Modulation of the preamble, message and pilot
[Qsymb_pre] = modulation(params,Preamble,params.modulation);      % Preamble modulation
[Qsymb_data] = modulation(params,bits_data,params.modulation);    % Message modulation
[Qsymb_pilot] = modulation(params,bits_pilot,params.modulation);  % Pilot modulation             

% 3. OFDM Transmitter: 
[signal_tx] = transmitter_Test(params, Qsymb_pre, Qsymb_data, Qsymb_pilot);

preambleLCP = signal_tx(:,1:params.Q+params.LCP);

progress_indx = 0;

for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(params.SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = params.SNR_list(SNR_idx);

        % 4. Channel propagation: 
        signal_rx = channel_propagation_test(params,signal_tx,SNR,STO,CFO,Nr);
        
        % 5. CFO and STO estimation:
        [STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx,1/(params.nPreamble+params.nData+params.N_zeros));
 
        % Average over the antennas
        STO_estimated = round(mean(STO_estimated,'all'));
        CFO_estimated = mean(CFO_estimated,'all');
        
        % STO correction
        signal_rx = signal_rx(:,STO_estimated+ones(size(signal_rx,1),1):STO_estimated+Nsymb*ones(size(signal_rx,1),1));
        
        % CFO correction
        T = 1/params.B;
        n = 1:1:Nsymb;
        phi = exp(1i*CFO_estimated*T*n);    
        signal_rx = signal_rx.*phi;

        preamble = Qsymb_pre(1:params.nActiveQ);

        % 6. OFDM Receiver:
        [hz,Qsymb_rx] = receiver_Test(params,signal_rx, preamble,Qsymb_pilot);

        % 7. Demodulation:
        bits_rx = demodulation(params,Qsymb_rx,params.modulation);
        
        % compute BER
        bits_tx = bits_data;
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
ber_theo = berawgn(params.SNR_list,'p�&am',2^(params.Nbps));
semilogy(params.SNR_list,ber_theo,'--');
legend('simulated','theoretical');
xlabel('SNR dB');ylabel('Probability of error');
xlim([-5 10]);
s = sprintf('SISO OFDM modulation BER in presence of AWGN noise \n QPSK case');
title(s);