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

function symb_rx = receiver(params,signal_rx,Nsymb_ofdm, Preamble)

   
    

    % S/P conversion
    s = reshape(signal_rx,length(signal_rx)/Nsymb_ofdm,Nsymb_ofdm);
    
    % CP removal
    s=s(params.ofdm.cp_L+1:end,:);
    
    %FFT
    S=fft(s(:,1:end),params.ofdm.N_subcrr);
    
    
    % Inactive subcarriers removal
    S = S((params.ofdm.N_inactive_subcrr-1)/2 + 1:end - (params.ofdm.N_inactive_subcrr-1)/2 ,:);
    N_active_subcrr = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;
    S = vertcat(S(1:(N_active_subcrr-1)/2,:),S(end - (N_active_subcrr-1)/2:end,:));  
    
    %Channel estimation
    lambda=diag(Preamble(:,2));    
    H= lambda\S(:,2);  
    h=ifft(H);
    h(1,257:end)=0;
    H=fft(h);
    %Channel equalization
    Hcirc = toeplitz(H, [H(1,1); zeros(N_active_subcrr-1,1)]);
    Hm = ones(N_active_subcrr, Nsymb_ofdm).*H;
%     Hm=zeros(N_active_subcrr, Nsymb_ofdm);
%     for i=1:1:Nsymb_ofdm
%     Hm(:,i)=H;
%     end
    S=S./Hm;
  
    
  
    % P/S conversion
    %symb_rx = reshape(S,params.ofdm.N_subcrr*Nsymb_ofdm,1);
    symb_rx = reshape(S,[],1);
    
    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end


