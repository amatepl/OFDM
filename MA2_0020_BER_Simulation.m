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
params.Nbps = 2;                        % Modulation order
params.modulation = 'qpsk';
params.Q = 2048;
params.nActiveQ = params.Q-410;
params.ActiveQIndex = [2:params.nActiveQ/2 params.Q-params.nActiveQ/2:params.Q];
params.nData = 30;
NsimPerSNR = 10;                        % number of simulations per SNR value

%% --- Local parameters
STO = 0;                                % Time offset (switching unit vector)
% delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
CFO = 0;           % Carrier frequency offset in ppm
%Nr = 1;                                 % number of receivers

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

BER_i=zeros(NsimPerSNR,length(params.SNR_list),4);

progress_indx = 0;
MSTO = zeros(NsimPerSNR,length(params.SNR_list));
MCFO = zeros(NsimPerSNR,length(params.SNR_list));
for i=1:4
for sim_idx = 1:NsimPerSNR
    for SNR_idx = 1:length(params.SNR_list)
        progress_indx =  progress_indx + 1;
        SNR = params.SNR_list(SNR_idx);

        % 4. Channel propagation: 
        signal_rx = channel_propagation_test(params,signal_tx,SNR,STO,CFO,i);
        
        % 5. CFO and STO estimation:
        [STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx,1);
        % For SISO only
        %MSTO(sim_idx,SNR_idx) = immse(STO,STO_estimated);
        %MCFO(sim_idx,SNR_idx) = immse(CFO*params.Fc,CFO_estimated);
        % Average over the antennas
%         STO_estimated = round(mean(STO_estimated,'all'));
%         CFO_estimated = mean(CFO_estimated,'all');
         
%         % STO correction
%         signal_rx = signal_rx(:,STO_estimated+ones(size(signal_rx,1),1):STO_estimated+Nsymb*ones(size(signal_rx,1),1));
%          
%         % CFO correction
%         T = 1/params.B;
%         n = 1:1:Nsymb;
%         phi = exp(-1i*CFO_estimated*T*n);    
%         signal_rx = signal_rx.*phi;

        preamble = Qsymb_pre(1:params.nActiveQ);
        
        % 6. OFDM Receiver:
        [hz,Qsymb_rx] = receiver_Test(params,signal_rx, preamble,Qsymb_pilot);  
        
        % 7. Demodulation:
        if params.N_pilots > 0
            bits_rx = demodulation(params,Qsymb_rx,params.modulation);
        else
            bits_rx = demodulation(params,Qsymb_rx(2*params.nActiveQ+1:end),params.modulation);
        end
        
        % compute BER
        bits_tx = bits_data;
        bitErrorRate = sum(abs(bits_tx - bits_rx),'all');
        BER_i(sim_idx,SNR_idx,i) = bitErrorRate / length(bits_tx);
    end
    progress = round(progress_indx/NsimPerSNR/length(params.SNR_list)*100);
    disp([' Simulation executed @ ', num2str(progress), '%']);
    disp('');
end
end

% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------

disp('$$ Displaying results:');
BER1=BER_i(:,:,1);
BER2=BER_i(:,:,2);
BER3=BER_i(:,:,3);
BER4=BER_i(:,:,4);


figure;
semilogy(params.SNR_list,mean(BER1,1),'-');
hold on;
semilogy(params.SNR_list,mean(BER2,1),'-');
hold on;
semilogy(params.SNR_list,mean(BER3,1),'-');
hold on;
semilogy(params.SNR_list,mean(BER4,1),'-');
grid on; hold on;
% ber theoretical
ber_theo = berawgn(params.SNR_list,'psk',2^(params.Nbps),'non-diff');
semilogy(params.SNR_list,ber_theo,'--');
legend('Nr=1','Nr=2','Nr=3','Nr=4','Theoretical','location','best');
xlabel('SNR dB');ylabel('Probability of error');
xlim([-5 13]);
s = sprintf("SISO OFDM modulation BER Diversity gain with 5 MPC's");
title(s);

% BER_i_5 = load("BERCFO_5.mat").BER_i;
% BER_i_15 = load("BERCFO_15.mat").BER_i;
% BER_i_30 = load("BERCFO_30.mat").BER_i;
% figure;
% semilogy(params.SNR_list,mean(BER_i_5,1),'--o');
% hold on;
% semilogy(params.SNR_list,mean(BER_i_15,1),'--+');
% hold on;
% semilogy(params.SNR_list,mean(BER_i_30,1),'--s');
% hold on;
% semilogy(params.SNR_list,mean(BER_i,1),'--^');
% grid on; hold on;
% % ber theoretical
% ber_theo = berawgn(params.SNR_list,'psk',2^(params.Nbps),'non-diff');
% semilogy(params.SNR_list,ber_theo,'--');
% legend('N = 5','N = 15','N = 30','N = 60','Theoretical','location','best');
% xlabel('SNR dB');ylabel('Probability of error');
% xlim([-5 13]);
% s = sprintf("SISO OFDM modulation BER for different sizes of ofdm symbol");
% title(s);

% figure;
% plot(params.SNR_list,mean(MSTO,1));
% grid on;
% ylabel("MSE ToA");
% xlabel('SNR dB');
% title("ToA estimate Mean Square Error function of the SNR");

% figure;
% plot(params.SNR_list,mean(MCFO,1));
% grid on;
% ylabel("MSE CFO");
% xlabel('SNR dB');
% title("CFO estimate Mean Square Error function of the SNR");