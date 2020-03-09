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

signal_rx = load('sig_rx.mat'); % load singal_rx
signal_rx = signal_rx.sig_rx;


% -------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------


% 4. OFDM Receiver:
[STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx);

%Average over the antennas
STO_estimated = round(mean(STO_estimated,'all'));
CFO_estimated = mean(CFO_estimated,'all');

T = 1/params.B;
n = 1:1:size(signal_rx,2);
phi = exp(1i*CFO_estimated*T*n);    
signal_rx = signal_rx.*phi;

% signal_rx = horzcat(signal_rx(STO_estimated+1:end),zeros(1,STO_estimated));

Nsymb = (params.Q+params.LCP)*(params.nData+params.nPreamble);
signal_rx = signal_rx(STO_estimated+1:STO_estimated+Nsymb);

Qsymb_pre = params.PreambleSymbols;
Qsymb_pre = Qsymb_pre(params.ActiveQIndex); % Removing inactive subcarriers

Qsymb_pilot = 0;

Qsymb_rx = receiver_Test(params,signal_rx,params.nData, Qsymb_pre,Qsymb_pilot);

% 5. Demodulation:
bits_rx = demodulation(params,Qsymb_rx,'bpsk');


% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------
% 
% bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);
% disp('$$ Displaying results:');
% disp(['BER:', num2str(bitErrorRate)]);
% 
% figure;
% subplot(1,2,1); plot(real(Qsymb_tx),imag(Qsymb_tx),'rx'); 
% title('Tx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
% subplot(1,2,2); plot(real(Qsymb_rx),imag(Qsymb_rx),'.'); 
% title('Rx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])

figure;
plot(real(Qsymb_rx),imag(Qsymb_rx),'.'); 
% title('Rx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])