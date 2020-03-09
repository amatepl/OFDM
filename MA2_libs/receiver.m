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

function symb_rx = receiver(params,signal_rx,Nsymb_ofdm, Preamble,pilot)

k=size(signal_rx,1);   
Hcomb = zeros(params.ofdm.N_active_subcrr,1);
for i=1:k  
    % S/P conversion
    %s = reshape(signal_rx,length(signal_rx)/Nsymb_ofdm,Nsymb_ofdm+2);
    s = reshape(signal_rx(i,:).',[],Nsymb_ofdm+2);
    
    % CP removal
    s=s(params.ofdm.cp_L+1:end,:);
    
    %FFT
    S=fft(s(:,1:end),params.ofdm.N_subcrr);
    
    
    % Inactive subcarriers removal
    S = S((params.ofdm.N_inactive_subcrr)/2 :end - (params.ofdm.N_inactive_subcrr)/2,:);
    N_active_subcrr = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;

    S = vertcat(S(1:(N_active_subcrr)/2,:),S(end - (N_active_subcrr)/2 +1:end,:));  
    
    %Channel estimation i frequency damain
%     lambda=diag(Preamble(:,2));    
    lambda=diag(Preamble(1:end - size(Preamble,1)/2,1));   

    H= lambda\S(:,2);  
    h=ifft(H); % refinement in time domain of f domain estimation
    h(1,257:end)=0;
    H=fft(h);
    %Channel equalization
    Hm = ones(N_active_subcrr, Nsymb_ofdm +2).*H;
    %S=S./Hm;

   
    %Channel estimation in time domain
    a=lambda'*S(:,2);
    ht = ifft(a,N_active_subcrr);
    Ht=fft(ht,N_active_subcrr); 
    %Channel equalization
    Htcirc = toeplitz(Ht, [Ht(1,1); zeros(N_active_subcrr-1,1)]);
    Htm = ones(N_active_subcrr, Nsymb_ofdm +2).*Ht;

    S=S.';%reshape(S,Nsymb_ofdm+2,N_active_subcrr);
    S=S*diag(Ht)';
    %S=reshape(S,N_active_subcrr,Nsymb_ofdm+2);
    S=S.';
  
    % CFO tracking
    S_pilots = S(:,3:end);
    
    S_pilots_1 = S_pilots(1:size(S_pilots,1)/2 ,:);
    S_pilots_2 = S_pilots(size(S_pilots,1)/2 +1:end,:);
    
    Hcomb_1 = Ht(1:size(Ht,1)/2 ,:);
    Hcomb_2 = Ht(size(Ht,1)/2 +1:end,:);
    
    S_pilots_1 = reshape(S_pilots_1,[],Nsymb_ofdm,params.ofdm.N_pilots/2);
    S_pilots_2 = reshape(S_pilots_2,[],Nsymb_ofdm,params.ofdm.N_pilots/2);
    
    Hcomb_1 = reshape(Hcomb_1,[],params.ofdm.N_pilots/2);
    Hcomb_2 = reshape(Hcomb_2,[],params.ofdm.N_pilots/2);
    
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
    
    Hcomb = vertcat(Hcomb_1,Hcomb_2);
    
    
    % Find impulse responses corresponding to pilot frequencies
    H_pilots_1 = Ht(1:size(Ht,1)/2 ,:);
    H_pilots_2 = Ht(size(Ht,1)/2 +1:end,:);
    
    H_pilots_1 = reshape(H_pilots_1,[],params.ofdm.N_pilots/2);
    H_pilots_2 = reshape(H_pilots_2,[],params.ofdm.N_pilots/2);
    
    H_pilots_1 = H_pilots_1(1,:);
    H_pilots_2 = H_pilots_2(end,:);
    
    H_pilots_1 = reshape(H_pilots_1,[],1);
    H_pilots_2 = reshape(H_pilots_2,[],1);
    
    H_pilots = vertcat(H_pilots_1,H_pilots_2);
    
    phi = conj(pilots_rx).*(H_pilots.*ones(size(pilots_rx)).*pilot);
    
    phi = angle(sum(phi,1));

  
    % P/S conversion
    %symb_rx = reshape(S,params.ofdm.N_subcrr*Nsymb_ofdm,1);
%     symb_rx = reshape(S,[],1);

    S_pilots = S_pilots.*(kron(exp(1i*phi),ones(size(S_pilots,1),1))); 
    
    Ssum =+ S_pilots;
    Hsum =+ abs(Hcomb).^2;

%     S(:,1:2) = S(:,1:2).*(kron(exp(1i*phi),ones(size(S,1),1)));
    
end    
    Scomb= Ssum./(Hsum.*ones(N_active_subcrr-params.ofdm.N_pilots, Nsymb_ofdm));
    Scomb = reshape(Scomb,[],1);
    
    %S_pilots = reshape(S_pilots,[],1);
    
%     S = reshape(S(:,1:2),[],1);     % Pilot symbols
    
  %  symb_rx = vertcat(S,S_pilots).*exp(1i*phi);
%     symb_rx = vertcat(S,S_pilots);
    symb_rx = Scomb;

    % ---------------------------------------------------------------------
    % 'simple_ofdm_Tx': Implements a simple ofdm transmitter: S/P, IFFT, CP
    % P/S. It doesn't consider active/inactive subcarriers.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % symb_rx = simple_ofdm_Rx(params,signal_rx,Nsymb_ofdm);
end