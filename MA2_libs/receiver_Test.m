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

function symb_rx = receiver_Test(params,signal_rx,Nsymb_ofdm, Preamble,pilot)
    k=size(signal_rx,1);
%     k = 1;

    figure, hold on;
    grid on;
    title("signal rx with CFO and STO correction applied");
    plot(abs(signal_rx));


    for i=1:k
        signalrx = signal_rx(i,1:end-mod(size(signal_rx,2),32));

        % S/P conversion
        s = reshape(signalrx,[],Nsymb_ofdm+params.nPreamble);

        % CP removal
        s=s(params.LCP+1:end,:);

        %FFT
        S=fft(s(:,1:end),params.Q);

%         fq = -params.Q/2:1:params.Q/2-1;
%         plot(fq,abs(ifftshift(S(:,1))));

        % Inactive subcarriers removal

        S = S(params.ActiveQIndex,:);
        

%         figure, hold on;
%         fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%         plot(fq,abs(fftshift(S(:,1))));
%         grid on;
%         title("signal rx after inactive subcarriers removal");

        %Channel estimation i frequency damain
        lambda=diag(Preamble);    
    %     lambda=diag(Preamble(1:end - size(Preamble,1)/2,1));   

        H= lambda\S(:,2);  
        h=ifft(H);
        h(1,params.LCP:end)=0;
        H=fft(h);
        H = H/sqrt(params.nActiveQ);
        %Channel equalization
        Hcirc = toeplitz(H, [H(1,1); zeros(params.nActiveQ-1,1)]);
        Hm = ones(params.nActiveQ, params.nData +2).*H;
        Sm = S./Hm;

%         figure, hold on;
%         fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%         plot(fq,abs(fftshift(H)));
%         grid on;
%         title("Frequency domain channel estimation");

        %Channel estimation in time domain
        a=lambda'*S(:,2);

        ht = dftmtx(params.nActiveQ)'*(lambda')*lambda*dftmtx(params.nActiveQ);

        ht = ht\ifft(a,params.nActiveQ);
        Ht=fft(ht(1:end,1),params.nActiveQ); 
        Ht = Ht/sqrt(params.nActiveQ);

        %Channel equalization

%         figure, hold on;
%         fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%         plot(fq,abs(fftshift(Ht)));
%         grid on;
%         title("Time domain channel estimation");
        
%         figure, hold on;
%         fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
%         ht = ifftshift(ifft(Ht,params.nActiveQ));
%         plot(fq,10*log10(ht));
%         grid on;
%         title("Time domain channel estimation (Time domain)");

        N_pilots = 126;

        Htcirc = toeplitz(Ht, [Ht(1,1); zeros(params.nActiveQ-1,1)]);
        Htm = ones(params.nActiveQ, params.nData +2).*Ht;

        S=S.';  % reshape(S,Nsymb_ofdm+2,N_active_subcrr);
        S=S*diag(Ht)';
        %S=reshape(S,N_active_subcrr,Nsymb_ofdm+2);
        S=S.';
        
        % CFO tracking
         S_pilots = S(:,3:end);
        
        S_pilots_1 = S_pilots(1:size(S_pilots,1)/2 ,:);
        S_pilots_2 = S_pilots(size(S_pilots,1)/2 +1:end,:);
        
        Hcomb_1 = Ht(1:size(Ht,1)/2 ,:);
        Hcomb_2 = Ht(size(Ht,1)/2 +1:end,:);
        
        S_pilots_1 = reshape(S_pilots_1,[],Nsymb_ofdm,N_pilots/2);
        S_pilots_2 = reshape(S_pilots_2,[],Nsymb_ofdm,N_pilots/2);
              
        Hcomb_1 = reshape(Hcomb_1,[],N_pilots/2);
        Hcomb_2 = reshape(Hcomb_2,[],N_pilots/2);
        
        % Extracting the recieved pilots
        pilots_rx_1 = S_pilots_1(1,:,:);
        pilots_rx_2 = S_pilots_2(end,:,:);
        
        pilots_rx_1 = reshape(pilots_rx_1,[],Nsymb_ofdm);
        pilots_rx_2 = reshape(pilots_rx_2,[],Nsymb_ofdm);
        
        pilots_rx = vertcat(pilots_rx_1,pilots_rx_2);   % Received pilots
        
        S_pilots_1 = S_pilots_1(2:end,:,:);
        S_pilots_2 = S_pilots_2(1:end-1,:,:);
        
        Hcomb_1 = Hcomb_1(2:end,:);
        Hcomb_2 = Hcomb_2(1:end-1,:);
        
        S_pilots_1 = reshape(S_pilots_1,[],Nsymb_ofdm);
        S_pilots_2 = reshape(S_pilots_2,[],Nsymb_ofdm);
        
        Hcomb_1 = reshape(Hcomb_1,[],1);
        Hcomb_2 = reshape(Hcomb_2,[],1);
        
        S_pilots = vertcat(S_pilots_1,S_pilots_2);  % Message without pilots
        
        Hcomb = vertcat(Hcomb_1,Hcomb_2);   % Impulse response for message symbols subcarriers
        
        % Find impulse responses corresponding to pilot frequencies
%         H_pilots_1 = Ht(1:size(Ht,1)/2 ,:);
%         H_pilots_2 = Ht(size(Ht,1)/2 +1:end,:);
%         
%         H_pilots_1 = reshape(H_pilots_1,[],N_pilots/2);
%         H_pilots_2 = reshape(H_pilots_2,[],N_pilots/2);
%         
%         H_pilots_1 = H_pilots_1(1,:);
%         H_pilots_2 = H_pilots_2(end,:);
%         
%         H_pilots_1 = reshape(H_pilots_1,[],1);
%         H_pilots_2 = reshape(H_pilots_2,[],1);
%         
%         H_pilots = vertcat(H_pilots_1,H_pilots_2);
        
        phi = conj(pilots_rx)*pilot;
        phi = angle(sum(phi,1));
    
      
        % P/S conversion
    
        S_pilots = S_pilots.*(kron(exp(1i*phi),ones(size(S_pilots,1),1))); 
        
        Ssum =+ S_pilots;
        Hsum =+ abs(Hcomb).^2;
    end
    
    Scomb= Ssum./(Hsum.*ones(params.nActiveQ-N_pilots, Nsymb_ofdm));
    Scomb = reshape(Scomb,[],1);

    
    
    figure, hold on;
    fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
    plot(fq,abs(fftshift(H)));
    grid on;
    title("Frequency domain estiamtion");

    figure, hold on;
    plot(fq,abs(fftshift(Ht)));
    grid on;
    title("Time domain estiamtion")

    symb_rx = Scomb;


        % ---------------------------------------------------------------------
        % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
        % P/S. It doesn't consider active/inactive subcarriers.
        % IMPORTANT!!: Comment the next line when trying your implementation
        % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end


