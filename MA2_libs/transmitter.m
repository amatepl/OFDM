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

function [signal_tx, Preamble] = transmitter(params,symb_tx,Nsymb_ofdm)
    
    % Serial to parallel converter
%     symb_tx_parallel = reshape(symb_tx,params.ofdm.N_subcrr,Nsymb_ofdm);
    symb_tx_parallel = reshape(symb_tx,[],Nsymb_ofdm);
    
    N_active_subcrr = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;
    % Inactive subcarriers removal
%     symb_tx_parallel(1:(params.ofdm.N_inactive_subcrr-1)/2,:) = 0;
%     symb_tx_parallel(end - (params.ofdm.N_inactive_subcrr-1)/2 +1:end,:) = 0;
%     symb_tx_parallel(params.ofdm.N_subcrr/2,:) = 0;     % DC

    symb_tx_parallel = vertcat(symb_tx_parallel(1:(N_active_subcrr-1)/2,:),...
                                zeros(1,Nsymb_ofdm),...
                                symb_tx_parallel(end - (N_active_subcrr-1)/2:end,:));

    symb_tx_parallel = vertcat(zeros((params.ofdm.N_inactive_subcrr-1)/2,Nsymb_ofdm),...
                               symb_tx_parallel,...
                               zeros((params.ofdm.N_inactive_subcrr-1)/2,Nsymb_ofdm));

    
    Preamble = symb_tx_parallel(:,1:2);
    
    
    % IFFT
    symb_tx_parallel = ifft(symb_tx_parallel,[],1);
   % symb_tx_parallel = ifftshift(symb_tx_parallel);
    
    % Cyclic prefix addition
    symb_tx_parallel = vertcat(symb_tx_parallel(end-params.ofdm.cp_L+1:end,:),symb_tx_parallel);
    
    % Parallel to serial converter
    signal_tx = reshape(symb_tx_parallel,1,[]);
    
    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % signal_tx = simple_ofdm_Tx(params,symb_tx,Nsymb_ofdm);
end
