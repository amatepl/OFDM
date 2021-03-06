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


function [Hsave,symb_rx] = receiver_Test(params,signal_rx,Preamble,pilot)
    Nsymb_ofdm = params.nData;
    k=size(signal_rx,1); 
    fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
    Hsave = zeros(params.nActiveQ,size(signal_rx,1));
   
    % S/P conversion
    s = reshape(signal_rx.',[],Nsymb_ofdm+params.nPreamble,size(signal_rx,1));
    
    % CP removal
    s = s(params.LCP+1:end,:,:);
    
    %FFT
    S=fft(s,params.Q,1);
    
    %fq = -params.Q/2:1:params.Q/2-1;
    %plot(fq,abs(ifftshift(S(:,1))));
    
    % Inactive subcarriers removal
    S = S(params.ActiveQIndex,:,:);
    
    %Channel estimation i frequency damain
    lambda=diag(Preamble);    
    %lambda=diag(Preamble(1:end - size(Preamble,1)/2,1));   

    H= lambda\S(:,2);  
    h=ifft(H); % refinement in time domain of f domain estimation
    h(1,params.LCP+1:end)=0;
    H=fft(h);
    %Channel equalization
    Hm = ones(params.nActiveQ, Nsymb_ofdm +2).*H;
    %S=S./Hm;

   
    %Channel estimation in time domain
    a=conj(Preamble).*S(:,2,:);
    ht = ifft(a,params.nActiveQ,1);
    Ht=fft(ht,params.nActiveQ,1); 
    Hsave = reshape(Ht,params.nActiveQ,size(Ht,3));

%     plot(fq,abs(Ht));
%     grid on;
    
    %Channel equalization
    Htcirc = toeplitz(Ht, [Ht(1,1); zeros(params.nActiveQ-1,1)]);
    Htm = ones(params.nActiveQ, Nsymb_ofdm +2).*Ht;

    % match filter
    S=S.*conj(Ht);
  
    % CFO tracking
    S_pilots = S(:,3:end,:);
    if params.N_pilots > 0

        for i = 1:size(S_pilots,3)
            S_inter = S_pilots(:,:,i);
            S_pilots_1 = S_inter(1:size(S_inter,1)/2 ,:);
            S_pilots_2 = S_inter(size(S_inter,1)/2 +1:end,:);

            Hcomb_1 = Ht(1:size(Ht,1)/2 ,:,i);
            Hcomb_2 = Ht(size(Ht,1)/2 +1:end,:,i);

            S_pilots_1 = reshape(S_pilots_1,[],Nsymb_ofdm,params.N_pilots/2);
            S_pilots_2 = reshape(S_pilots_2,[],Nsymb_ofdm,params.N_pilots/2);

            Hcomb_1 = reshape(Hcomb_1,[],params.N_pilots/2);
            Hcomb_2 = reshape(Hcomb_2,[],params.N_pilots/2);

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

            S_inter = vertcat(S_pilots_1,S_pilots_2);  % Message without pilots

            Hcomb = vertcat(Hcomb_1,Hcomb_2);   % Impulse response for message symbols subcarriers


            % Find impulse responses corresponding to pilot frequencies
            H_pilots_1 = Ht(1:size(Ht,1)/2 ,:,i);
            H_pilots_2 = Ht(size(Ht,1)/2 +1:end,:,i);

            H_pilots_1 = reshape(H_pilots_1,[],params.N_pilots/2);
            H_pilots_2 = reshape(H_pilots_2,[],params.N_pilots/2);

            H_pilots_1 = H_pilots_1(1,:);
            H_pilots_2 = H_pilots_2(end,:);

            H_pilots_1 = reshape(H_pilots_1,[],1);
            H_pilots_2 = reshape(H_pilots_2,[],1);

            H_pilots = vertcat(H_pilots_1,H_pilots_2);  % Impulse response of pilot's subcarriers 

            phi = conj(pilots_rx).*(H_pilots.*ones(size(pilots_rx)).*pilot);

            phi = angle(sum(phi,1));

            S_inter = S_inter.*(kron(exp(1i*phi),ones(size(S_inter,1),1))); 

            Hsum =+ abs(Hcomb).^2;
            Ssum =+ S_inter;
        end
        Scomb= Ssum./(Hsum.*ones(params.nActiveQ-params.N_pilots, Nsymb_ofdm));
    else 
        Ssum = sum(S,3);
        Hsum = sum(abs(Ht).^2,3);
        Scomb= Ssum./(Hsum.*ones(params.nActiveQ, Nsymb_ofdm+params.nPreamble));
    end
    
    % P/S conversion
    Scomb = reshape(Scomb,[],1);
    
    %     figure, hold on;
    %     fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;
    %     plot(fq,abs(fftshift(H)));
    %     grid on;
    %     title("Frequency domain estiamtion");

    %     figure, hold on;
    %     plot(abs(Hsum));
    %     grid on;
    %     title("Time domain estiamtion")

    %     figure, hold on;
    %     plot(abs(ifft(Hsum)));
    %     grid on;
    %     title("Time domain estiamtion - h(t)")
    
    symb_rx = Scomb;

    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end
