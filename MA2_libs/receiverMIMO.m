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

function [symb_rx1, symb_rx2] = receiverMIMO(params,signal_rx,W,H1,H2,Nr)

% fq = -params.nActiveQ/2:1:params.nActiveQ/2-1;

% figure, hold on;
% grid on;
% title("Time domain estimation")

    S = signal_rx(:,:,1:Nr);
    
    % Inactive subcarriers removal
    S = S(params.ActiveQIndex,:,:);
    
    S = permute(S,[3 2 1]);
    
    W = permute(W,[2 1 3]);
    
    symb_rx1 = zeros(params.nActiveQ,size(S,2)); 
    symb_rx2 = zeros(params.nActiveQ,size(S,2));
    
    % Zero forcing equalizer
    
    % The following loop correcponds to W*Y for every symbol
    for i = 1:size(S,2)
        symb_rx1(:,i) = sum(W(:,1,:).*S(:,i,:),1);
        symb_rx2(:,i) = sum(W(:,2,:).*S(:,i,:),1);
    end
    
%     figure, hold on;
%     plot(abs(symb_rx1(:,1)));
%     title('Symb rx1 eqalized');
%     
%     figure, hold on;
%     plot(abs(symb_rx2(:,1)));
%     title('Symb rx2 eqalized')
    
%     keyboard
    
     % Zero forcing equalizer
%     H1 = fftshift(H1,2);
%     H1 = permute(H1,[2 3 1]);
%     H1 = H1(params.ActiveQIndex,:,:);
%     
%     H2 = fftshift(H2,2);
%     H2 = permute(H2,[2 3 1]);
%     H2 = H2(params.ActiveQIndex,:,:);
%      
%     S1 = S./H1;
%     S2 = S./H2;
% keyboard
%     % match filter
% %     S=S.*conj(Ht);
%         
%     Ssum1 = sum(S1,3);
%     Ssum2 = sum(S2,3);
%     Hsum = sum(abs(Ht).^2,3);

%     Scomb= Ssum./(Hsum.*ones(params.nActiveQ, Nsymb_ofdm+params.nPreamble));
%     Scomb1= Ssum1;
%     Scomb2= Ssum2;
    
    % P/S conversion
%     Scomb1 = reshape(Scomb1,[],1);
%     Scomb2 = reshape(Scomb2,[],1);
    
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
    
%     symb_rx1 = Scomb1;
%     symb_rx2 = Scomb2;

    symb_rx1 = reshape(symb_rx1,[],1);
    symb_rx2 = reshape(symb_rx2,[],1);
end