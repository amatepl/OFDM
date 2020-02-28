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

clear; close all; clc;
addpath(genpath('MA2_libs'));           % add libraries

cfg = load('MA2_lab_parameters.mat');   % load configFile
params = cfg.params;                    % get the set of parameters
dispConfigFile(params);                 % display the parameters

% --- Local parameters
SNR = 20;           % SNR in dB
STO = 30;
% Nsymb_ofdm = 2;     % number OFDM symbols to transmit
Nsymb_ofdm = params.ofdm.data_L;     % number OFDM symbols to transmit
Nbits = Nsymb_ofdm * (params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr) * params.modulation.Nbps;



% -------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------

% 1. QAM Modulation.
[bits_tx,Qsymb_tx, Preamble] = modulation(params,Nbits);

% 2. OFDM Transmitter: 
[signal_tx, Preamble_mod] = transmitter(params,Qsymb_tx,Nsymb_ofdm);

% 3. Channel propagation: 
signal_rx = channel_propagation(params,signal_tx,SNR,STO);

[STO_estimated, ~] = estimationSTOCFO(params,signal_rx);

% 4. OFDM Receiver:
Qsymb_rx = receiver(params,signal_rx,Nsymb_ofdm, Preamble_mod);

% 5. Demodulation:
bits_rx = demodulation(params,Qsymb_rx);


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