% Opera - ULB
% February 2020
%
% ToDo: 
%       - Multipath Channel propagation.
%       - STO
%       - AWGN
%
% Inputs: 
%       params      : MA2 parameters. 
%       signal_tx   : transmitted signal in time domain. dim = (1,Nsamples)
%       SNR         : desired SNR in dB. dim = (1,1)
% Outputs:
%       signal_rx   : received signal at the output of the channel. 
%                     dim = (1,Nsamples)

function signal_rx = channel_propagation4(params,signal_tx,H,SNR,Nr)

    H = fftshift(H,2);
%     h = ifft(H,[],2);
   
    H = permute(H,[2 3 1]); % Tranforms H into a 2048x1x4 matrix
    
    signal_rx = signal_tx.*H;   % Multiplication in frequency domain of the channel and the signal
    
     % 2. Noise computation
    transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
    noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
    noise_var = noise_energy/(length(signal_tx(:)));   % variance of noise to be added
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
%     noise = noise_std*(randn(length(signal_tx),Nr)+1i*randn(length(signal_tx),Nr));      % noise
    noise = noise_std*(randn(size(signal_tx))+1i*randn(size(signal_tx)));      % noise

    
    
 
    
%     figure, hold on;
%     plot(abs(H(1,:)));
    
%     h = ifft((H(:,params.ActiveQIndex)),[],2);
%     figure, hold on;
%     plot(abs(h(1,:).'));
    

%     signal_rx = signal_tx;
    
%     figure,hold on;
%     plot(abs(fft(h(1,:),2048)));
    
%     for i = 1:Nr
%         signal_rx(i,:) = conv(signal_tx,h(i,:),'same');
% %         signal_rx(i,:) = conv(signal_tx,1,'same');
%     end

%     figure, hold on;
%     plot(abs((H(1,:)).'));


%     h = h(:,1:params.LCP);
%     signal_rx = conv2(1,signal_tx,h);
%     N = (params.nData + params.nPreamble)*(params.Q+params.LCP);
%     signal_rx = signal_rx(:,1:N);
    
   
    
    
   
    
    signal_rx = signal_rx+noise; 
    
end