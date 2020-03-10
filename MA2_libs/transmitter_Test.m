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

function [signal_tx] = transmitter_Test(params, symb_pre,symb_tx, symb_pilot)
    
    % Serial to parallel converter
%     symb_tx_parallel = reshape(symb_tx,params.ofdm.N_subcrr,Nsymb_ofdm);

    N_pilots = 126;

    symb_tx = reshape(symb_tx,[],params.nData);
    
    % Add pilots
    symb_tx_1 = symb_tx(1:size(symb_tx)/2,:);
    symb_tx_2 = symb_tx(size(symb_tx)/2+1:end,:);
    
    symb_tx_1 = reshape(symb_tx_1,[],params.nData,N_pilots/2);
    symb_tx_2 = reshape(symb_tx_2,[],params.nData,N_pilots/2);
    
    symb_tx_1 = padarray(symb_tx_1,[1 0 0],symb_pilot,'pre');
    symb_tx_2 = padarray(symb_tx_2,[1 0 0],symb_pilot,'post');
    
    symb_tx_1 = reshape(symb_tx_1,[],params.nData);
    symb_tx_2 = reshape(symb_tx_2,[],params.nData);
    
    symb_tx = vertcat(symb_tx_1,symb_tx_2);

    symb_pre = reshape(symb_pre,[],2);
    symb_tx_parallel = horzcat(symb_pre,symb_tx);
%     symb_tx_parallel = reshape(symb_tx,[],Nsymb_ofdm);
    
    
    
    
%     N_active_subcrr = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;
    N_inactive_subcrr = params.Q - params.nActiveQ;
    
    % Inactive subcarriers removal

% BasicVariableLocationPilot = {0, 39, 261,330,342,351,522,522,645,651,708,...
%                               726,756,792,849,855,918,21017,1143,1155,1158,...
%                               1185,1206,1260, 1407,1419,1428,1461,1530,1545,...
%                               1572,1701,1702};
% Pilot placement equation k = 3 * L + 12 * P 
% 
% k ∈ Indices from 0 to the number of Overall Usable Carriers minus 1
% L ∈ 0 .. 3 denotes the symbol number with a cyclic period of 4
% P ≥ 0 is an integer number
    % pilotQAM

%     continuesLocationPilots = padarray([1 0 0 0 0 0 0 0 0 0 0],[0 params.ofdm.N_active_subcrr - 11],'circ','post');
%     continuesLocationPilots = continuesLocationPilots * symb_pilot;
    
   
    

%     symb_tx_parallel = vertcat(symb_tx_parallel(1:(params.nActiveQ)/2,:),...
%                                 zeros(1,params.nData+params.nPreamble),...
%                                 symb_tx_parallel(end - (params.nActiveQ)/2 +1:end,:));
% 
%     symb_tx_parallel = vertcat(zeros((N_inactive_subcrr)/2 -1,params.nData+params.nPreamble),...
%                                symb_tx_parallel,...
%                                zeros((N_inactive_subcrr)/2,params.nData+params.nPreamble));

    
    inactSubRem = zeros(params.Q,params.nData + params.nPreamble);
    
    inactSubRem(params.ActiveQIndex,:) =  symb_tx_parallel;
    
    symb_tx_parallel = inactSubRem;
    
    % IFFT
    symb_tx_parallel = ifft(symb_tx_parallel,[],1);
   % symb_tx_parallel = ifftshift(symb_tx_parallel);
    
    % Cyclic prefix addition
    symb_tx_parallel = vertcat(symb_tx_parallel(end-params.LCP+1:end,:),symb_tx_parallel);
    
    % Parallel to serial converter
    signal_tx = reshape(symb_tx_parallel,1,[]);
    
    % Add zeros
    
    
    
    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % signal_tx = simple_ofdm_Tx(params,symb_tx,Nsymb_ofdm);
end
