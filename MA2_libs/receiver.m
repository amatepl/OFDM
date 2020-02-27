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
%     CP_length=(length(signal_rx)-params.ofdm.N_subcrr*Nsymb_ofdm)/Nsymb_ofdm;
%     s=s(CP_length+1:end,:);

    s=s(params.ofdm.cp_L+1:end,:);
    
    %FFT
    S=fft(s(:,1:2),params.ofdm.N_subcrr);
%     S = fftshift(S);

    figure, hold on;
    plot(abs(S));
    
    %Channel estimation
    H= Preamble(:,1)./Preamble(:,2);
    
    %Equalization: inversion
    S2=zeros(params.ofdm.N_subcrr,Nsymb_ofdm);
    for i=1:1:Nsymb_ofdm
        S2(:,i)=S(:,i)./H;
        for j=1:params.ofdm.N_subcrr
            if(isnan(S2(j,i)))
                S2(j,i)=0;
            end
        end   
    end
    
    % P/S conversion
    symb_rx = reshape(S,params.ofdm.N_subcrr*Nsymb_ofdm,1);
    
    
    
    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end


