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

cfg = load('MA2_lab_parameters.mat');   % load configFile
params = cfg.params;                    % get the set of parameters
dispConfigFile(params);                 % display the parameters

%% --- Local parameters
SNR = 20;                           % Wanted SNR in dB
STO = 100;                           % Time offset (switching unit vector)
% delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
CFO = 10e-6;                        % Carrier frequency offset
% Nsymb_ofdm = 2;                   % number OFDM symbols to transmit
Nsymb_ofdm = params.ofdm.data_L;    % number OFDM symbols to transmit
Nr = 2;                             % number of receivers

Nbps = params.modulation.Nbps;      % QAM modulation
%Nbps = 1;                          % BPSK modulation
% Number of bits knowing the inactive subcarriers and the number of pilots
Nbits = Nsymb_ofdm * (params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr- params.ofdm.N_pilots) * Nbps;


%% ------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------

% 1. Message, preamble and pilot construction
[Preamble, bits_data, bits_pilot] = build_message(params,Nbits,Nbps);

%bits_tx = vertcat(Preamble,bits_data);
 bits_tx = bits_data;
 
% 2. Modulation of the preamble, message and pilot
[Qsymb_pre] = modulation(Nbps,Preamble,'qpsk');      % Preamble modulation
[Qsymb_data] = modulation(Nbps,bits_data,'qpsk');    % Message modulation
[Qsymb_pilot] = modulation(Nbps,bits_pilot,'qpsk');  % Pilot modulation
% OFDM symbols [2 x preamble + message,1]
Qsymb_tx = vertcat(Qsymb_pre,Qsymb_data);             

% 3. OFDM Transmitter: 
[signal_tx] = transmitter(params, Qsymb_pre, Qsymb_data, Qsymb_pilot, Nsymb_ofdm);

% 4. Channel propagation: 
signal_rx = channel_propagation(params,signal_tx,SNR,STO,CFO,Nr);

% 4. OFDM Receiver:
[STO_estimated, CFO_estimated] = estimationSTOCFO(params,signal_rx);

%Average over the antennas
STO_estimated = round(mean(STO_estimated,'all'));
CFO_estimated = mean(CFO_estimated,'all');

%rough CFO correction
T = 1/params.ofdm.B;
n = 1:1:size(signal_rx,2);
phi = exp(1i*CFO_estimated*T*n);
signal_rx = signal_rx.*phi;

%rough STO correction
signal_rx = [signal_rx(:,STO_estimated+1:end),signal_rx(:,1:STO_estimated)];

Qsymb_rx = receiver(params,signal_rx,Nsymb_ofdm, Qsymb_pre,Qsymb_pilot);

% 5. Demodulation:
bits_rx = demodulation(params,Qsymb_rx,'qpsk');


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