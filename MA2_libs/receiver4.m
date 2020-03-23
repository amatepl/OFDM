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

function [hz,symb_rx] = receiver4(params,signal_rx,Nsymb_ofdm, Preamble)

k=size(signal_rx,3); 
% fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
hz = zeros(k,1);

% figure, hold on;
% grid on;
% title("Time domain estimation")

    S = signal_rx;
    
    % Inactive subcarriers removal
    S = S(params.ActiveQIndex,:,:);
    
    %Channel estimation i frequency damain
  
    H = S(:,2,:)./Preamble;  
    
    h=ifft(H,[],1); % refinement in time domain of f domain estimation
    h(params.LCP+1:end,1,:)=0;
    H=fft(h,[],1);
    %Channel equalization
    Hm = ones(params.nActiveQ, Nsymb_ofdm +2).*H;
   
    %Channel estimation in time domain
    a=conj(Preamble).*S(:,2,:);
    ht = ifft(a,params.nActiveQ,1);
    Ht=fft(ht,params.nActiveQ,1); 
%     hz(i) = Ht(400);
    
    %Channel equalization

    Htm = ones(params.nActiveQ, Nsymb_ofdm +2).*Ht;

    % match filter
%     S=S.*conj(Ht);
    
    % Zero forcing equalizer
    S=S./Ht;
        
    Ssum = sum(S,3);
    Hsum = sum(abs(Ht).^2,3);

    Scomb= Ssum./(Hsum.*ones(params.nActiveQ, Nsymb_ofdm+params.nPreamble));
    
    % P/S conversion
    Scomb = reshape(Scomb,[],1);
    
%     figure, hold on;
%     fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%     plot(fq,abs(fftshift(H(:,1,1))));
%     grid on;
%     title("Frequency domain estiamtion");
% 
%     figure, hold on;
%     fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%     plot(fq,abs(Ht(:,1,1)));
%     grid on;
%     title("Time domain estiamtion")
    
    symb_rx = Scomb;
end