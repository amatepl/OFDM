% Opera - ULB
% February 2020
%
% ToDo: 
%       - Implement the serial to parallel converter
%       - Consider Active/Inactive subcarriers.
%       - IFFT.
%       - Cyclic prefix Addition.
%       - Paralel to serial converter.
%
% Inputs: 
%       params      : MA2 parameters. 
%       symb_tx     : QAM symbols, generated in the modulation stage.
%                     dim = (N_qam_symb,1)
%       Nsymb_ofdm  : Number of OFDM symbols that are transmitted. 
%                     dim=(1,1)
% Outputs:
%       signal_tx   : signal in time domain to be transmitted. 
%                     dim=(1,Nsamples)
%
% Help: S/P example: 
%                       1 4 7
% 1 2 3 4 5 6 7 8 9 <=> 2 5 8
%                       3 6 9
%  dim(1,Nsymb_qam) <=> dim(N_subcrr,Nsymb_ofdm)

function [signal_tx] = transmitter4(params, symb_pre,symb_tx)

    % Serial to parallel converter
    symb_tx = reshape(symb_tx,[],params.nData);

    % Preamble addition
    symb_pre = reshape(symb_pre,[],2);
    symb_tx_parallel = horzcat(symb_pre,symb_tx);  
    
    % Adding inactive subcarriers
    inactSubRem = zeros(params.Q,params.nData + params.nPreamble);
    inactSubRem(params.ActiveQIndex,:) =  symb_tx_parallel;
    
    signal_tx = (inactSubRem);
    
end
