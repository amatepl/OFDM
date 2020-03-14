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

signal_tx = load('../ma2/sig_tx.mat').sig_tx; % load singal_rx

signal_rx_los = load('../ma2/ma2_g1_los_rx.mat'); % load singal_rx
signal_rx_los = signal_rx_los.ma2_g1_los_rx;

frame_size = (params.nPreamble+params.nData+params.N_zeros)*(params.Q+params.LCP);
Nsymb = (params.Q+params.LCP)*(params.nData+params.nPreamble);

%% ------------------------------------------------------------------------
% ------------------- OFDM Communication Chain ----------------------------
% -------------------------------------------------------------------------
Qsymb_pilot = -1;

preambleLCP = signal_tx(:,1).';
preamble = signal_tx(params.LCP+1:end,1:2);
Qsymb_pre=fft(preamble(:,1:end),params.Q);
Qsymb_pre = Qsymb_pre(params.ActiveQIndex,:);

data = signal_tx(params.LCP+1:end,3:2+params.nData);
Qsymb_data = fft(data(:,1:end),params.Q);
Qsymb_data = Qsymb_data(params.ActiveQIndex,:);
Qsymb_data_1 = Qsymb_data(1:size(Qsymb_data)/2,:);
Qsymb_data_2 = Qsymb_data(size(Qsymb_data)/2+1:end,:);
Qsymb_data_1 = reshape(Qsymb_data_1,[],params.nData,params.N_pilots/2);
Qsymb_data_2 = reshape(Qsymb_data_2,[],params.nData,params.N_pilots/2);
Qsymb_data_1 = Qsymb_data_1(2:end,:,:);
Qsymb_data_2 = Qsymb_data_2(1:end-1,:,:);
Qsymb_data_1 = reshape(Qsymb_data_1,[],params.nData);
Qsymb_data_2 = reshape(Qsymb_data_2,[],params.nData);
Qsymb_data = vertcat(Qsymb_data_1,Qsymb_data_2);
Qsymb_data = reshape(Qsymb_data,[],1);
bits_tx = demodulation(params,Qsymb_data,'bpsk');
Qsymb_preamble = reshape(Qsymb_pre,[],1);
Qsymb_tx = vertcat(Qsymb_preamble,Qsymb_data);

signal_rx = signal_rx_los(:,1:3*frame_size);

% 4. OFDM Receiver:
[STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx,frame_size,preambleLCP);

%Average over the antennas
STO_estimated = round(mean(STO_estimated,'all'));
CFO_estimated = mean(CFO_estimated,'all');

signal_rx = signal_rx(:,STO_estimated+1:STO_estimated+Nsymb);

T = 1/params.B;
n = 1:1:size(signal_rx,2);
phi = exp(1i*CFO_estimated*T*n);    
signal_rx = signal_rx.*phi;

[hz,Qsymb_rx] = receiver_Test(params,signal_rx,params.nData,Qsymb_pre(:,1).',Qsymb_pilot);

% 5. Demodulation:
bits_rx = demodulation(params,Qsymb_rx,'bpsk');


% -------------------------------------------------------------------------
% -------- Displaying results
% -------------------------------------------------------------------------
% 
bitErrorRate = sum(abs(bits_tx - bits_rx),'all') / length(bits_tx);
disp('$$ Displaying results:');
disp(['BER:', num2str(bitErrorRate)]);
% 
figure;
subplot(1,2,1); plot(real(Qsymb_tx),imag(Qsymb_tx),'rx'); 
title('Tx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])
subplot(1,2,2); plot(real(Qsymb_rx),imag(Qsymb_rx),'.'); 
title('Rx qam constellation');grid on; axis([-2,2,-2,2]);pbaspect([1 1 1])