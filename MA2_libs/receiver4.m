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
for i=1:k  
    
%     signalrx = signal_rx(i,1:end-mod(size(signal_rx,2),32));
%     
%     % S/P conversion
%     s = reshape(signalrx.',[],Nsymb_ofdm+params.nPreamble);
%     
%     % CP removal
%     s=s(params.LCP+1:end,:);
%     
%     %FFT
%     S=fft(s(:,1:end),params.Q);
    
    %fq = -params.Q/2:1:params.Q/2-1;
    %plot(fq,abs(ifftshift(S(:,1))));
    
%     figure, hold on;
%     plot(abs(S(:,2)));
%     title("S + inactQ");

    S = signal_rx(:,:,i);
    
    % Inactive subcarriers removal
    S = S(params.ActiveQIndex,:);
    
%     figure, hold on;
%     plot(abs(S(:,2)));
%     title("S");
    
%     figure, hold on;
%     plot((Preamble));
    
    %Channel estimation i frequency damain
    lambda=diag(Preamble);    

    H= lambda\S(:,2);  
    h=ifft(H); % refinement in time domain of f domain estimation
    h(1,params.LCP+1:end)=0;
    H=fft(h);
    %Channel equalization
    Hm = ones(params.nActiveQ, Nsymb_ofdm +2).*H;
    %S=S./Hm;

   
    %Channel estimation in time domain
    a=lambda'*S(:,2);
    ht = ifft(a,params.nActiveQ);
    Ht=fft(ht,params.nActiveQ); 
    hz(i) = Ht(400);
    
%     plot(fq,abs(Ht));
%     grid on;
    
    %Channel equalization
    Htcirc = toeplitz(Ht, [Ht(1,1); zeros(params.nActiveQ-1,1)]);
    Htm = ones(params.nActiveQ, Nsymb_ofdm +2).*Ht;

    S=S.';
    S=S*diag(Ht)';
    S=S.';
    
    Ssum =+ S;
    Hsum =+ abs(H).^2;
    
end    

    Scomb= Ssum./(Hsum.*ones(params.nActiveQ, Nsymb_ofdm+params.nPreamble));
    
    % P/S conversion
    Scomb = reshape(Scomb,[],1);
    
%     figure, hold on;
%     fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%     plot(fq,abs(fftshift(H)));
%     grid on;
%     title("Frequency domain estiamtion");

%     figure, hold on;
%     fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%     plot(fq,abs(Ht));
%     grid on;
%     title("Time domain estiamtion")
    
    symb_rx = Scomb;

    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end