% Opera - ULB
% February 2020
%
% ToDo: 
%       - STO, CFO estimation and correction.
%       - Implement S/P converter                V
%       - Cyclic prefix Removal.                 V
%       - FFT                                    V
%       - Channel estimation
%       - Channel Equalization
%       - CFO tracking
%       - Fine CFO compensation.
%       - P/S converter.
%
% Inputs: 
%       params      : MA2 parameters. 
%       signal_rx   : signal in time domain after the channel.
%                     dim = (1,Nsamples)
%       Nsymb_ofdm  : Number of OFDM symbols that are transmitted. 
%                     dim = (1,1)
% Outputs:
%       symb_rx     : QAM symbols. dim = (N_qam_symb,1)
%

function symb_rx = receiver(params,signal_rx,Nsymb_ofdm)
    % S/P conversion
    s = reshape(signal_rx,length(signal_rx)/Nsymb_ofdm,Nsymb_ofdm);
    % CP removal
    CP_length=(length(signal_rx)-2048*Nsymb_ofdm)/Nsymb_ofdm;
    s=s(CP_length+1:end,:);
    %FFT
    S=fft(s(:,1:end),2048);
    % P/S conversion
    symb_rx = reshape(S,2048*Nsymb_ofdm,1);
    
    
    
    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end


