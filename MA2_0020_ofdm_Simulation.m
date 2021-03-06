% Opera - ULB
% February 2020
%
% This script defines the main structure of an OFDM communication
% chain. The chain is divided into 5 stages:
%    1. Modulation  : Converts bits into quam symbols.
%    2. Transmitter : Converts quam symbols into OFDM frame.
%    3. Channel Propagation: Multipaths, STO, Noise, etc.
%    4. Receiver    : Converts OFDM fram into quam symbols, also implements
%                     estimation/correction of STO & CFO, channel 
%                     equalization, etc. 
%    5. Demodulation: Converts QAM symbols into bits.
%
% REMARK 1: Please ensure that the input/ouput of each stage matches to the
% format specified in the description given at the begining of each 
% function. Last but not least happy coding!
%
% Note: If you need to pass more parameters to the function you can add
% them to the variable: 'params' in the following way: 
%                  params.extra.newVarName = NewVarValue
% -------------------------------------------------------------------------

clear ; close all; clc;
addpath(genpath('MA2_libs'));           % add libraries

cfg = load('TestParam.mat');   % load configFile
params = cfg.TestParam;                    % get the set of parameters
dispConfigFile_Test(params);                 % display the parameters
params.N_pilots = 126;
params.N_zeros = 2;
params.Nbps = 2;                        % Modulation order
% params.Q = 2048;
% params.nActiveQ = params.Q-410;
% params.ActiveQIndex = [1:params.nActiveQ/2 params.Q-params.nActiveQ/2:params.Q-1];
params.modulation = 'qpsk';

%% --- Local parameters
SNR = 20;                           % Wanted SNR in dB
STO = 0;                           % Time offset (switching unit vector)
% delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
CFO = 30e-6;                        % Carrier frequency offset
% Nsymb_ofdm = 2;                   % number OFDM symbols to transmit
Nsymb_ofdm = params.nData;    % number OFDM symbols to transmit
Nr = 4;                             % number of receivers

% Nbps = params.modulation.Nbps;      % QAM modulation
% Number of bits knowing the inactive subcarriers and the number of pilots
Nbits = Nsymb_ofdm * (params.nActiveQ - params.N_pilots) * params.Nbps;
frame_size = (params.nPreamble+params.nData+params.N_zeros)*(params.Q+params.LCP);
Nsymb = (params.Q+params.LCP)*(params.nData+params.nPreamble);


%% ------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------

% 1. Message, preamble and pilot construction
[Preamble, bits_data, bits_pilot] = build_message_test(params,Nbits);

%bits_tx = vertcat(Preamble,bits_data);
bits_tx = bits_data;
 
% 2. Modulation of the preamble, message and pilot
[Qsymb_pre] = modulation(params,Preamble,params.modulation);      % Preamble modulation
[Qsymb_data] = modulation(params,bits_data,params.modulation);    % Message modulation
[Qsymb_pilot] = modulation(params,bits_pilot,params.modulation);  % Pilot modulation
% OFDM symbols [2 x preamble + message,1]
Qsymb_tx = vertcat(Qsymb_pre,Qsymb_data);             

% 3. OFDM Transmitter: 
[signal_tx] = transmitter_Test(params, Qsymb_pre, Qsymb_data, Qsymb_pilot);
preambleLCP = signal_tx(:,1:params.Q+params.LCP);

% % 4. Channel propagation: 
signal_rx = channel_propagation_test(params,signal_tx,SNR,STO,CFO,Nr);
% 4. OFDM Receiver:
[STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx,1);
%Average over the antennas
STO_estimated = round(mean(STO_estimated,'all'));
CFO_estimated = mean(CFO_estimated,'all');

signal_rx = signal_rx(:,STO_estimated+ones(size(signal_rx,1),1):STO_estimated+Nsymb*ones(size(signal_rx,1),1));

T = 1/params.B;
n = 1:1:Nsymb;
phi = exp(-1i*CFO_estimated*T*n);    
signal_rx = signal_rx.*phi;

% signal_rx = horzcat(signal_rx(STO_estimated+1:end),zeros(1,STO_estimated));
preamble = Qsymb_pre(1:params.nActiveQ);

% [hz,Qsymb_rx] = receiver_Test(params,signal_rx, preamble,Qsymb_pilot);
[Hsave,Qsymb_rx] = receiver_Test(params,signal_rx,preamble,Qsymb_pilot);

% 5. Demodulation:
if params.N_pilots > 0
    bits_rx = demodulation(params,Qsymb_rx,params.modulation);
else
    bits_rx = demodulation(params,Qsymb_rx(2*params.nActiveQ+1:end),params.modulation);
end


% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------

bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);
disp('$$ Displaying results:');
disp(['BER:', num2str(bitErrorRate)]);

figure;
subplot(1,2,1); plot(real(Qsymb_tx),imag(Qsymb_tx),'rx'); 
title('Tx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
subplot(1,2,2); plot(real(Qsymb_rx),imag(Qsymb_rx),'.'); 
title('Rx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);

% Qsymb_rx1 = load("MPC_const").Qsymb_rx;
% figure;
% subplot(1,2,1); plot(real(Qsymb_rx1),imag(Qsymb_rx1),'rx'); 
% title("QPSK constellation with 2 MPC's");grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
% subplot(1,2,2); plot(real(Qsymb_rx),imag(Qsymb_rx),'.'); 
% title('QPSK constellation with 2 simple equalizer');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
% bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);

% STX = fft(signal_rx(1:params.Q));
% STY = fft(signal_tx(1:params.Q));
% 
% figure;
% plot(abs(fftshift(STX)),'-o');
% hold on;
% plot(abs(fftshift(STY)),'-x');
% legend('TX side','RX side')
% grid on;
% title("Transmitted/Received OFDM symbol without noise and MPC's")