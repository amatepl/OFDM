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

function signal_rx = channel_propagationMIMO(params,signal_tx1,signal_tx2,H1,H2,SNR,Nr)

    H1 = fftshift(H1,2);
    H1 = permute(H1,[2 3 1]); % Tranforms H into a 2048x1x4 matrix
    
    signal_rx1 = signal_tx1.*H1;   % Multiplication in frequency domain of the channel and the signal
    
    H2 = fftshift(H2,2);
    H2 = permute(H2,[2 3 1]); % Tranforms H into a 2048x1x4 matrix
    
    signal_rx2 = signal_tx2.*H2;   % Multiplication in frequency domain of the channel and the signal
    
    signal_tx1 = signal_tx1+signal_tx2;
    
    signal_tx1 = ifft(signal_tx1,[],1);
    
     % 2. Noise computation
    transmitted_energy = (norm(signal_tx1(:)))^2;%.*vecnorm(H,2,1).^2;%*sum(vecnorm(H,2,1).^2,'all');           % energy of the signal
%     signal = reshape(signal_rx,params.Q*(params.nData + params.nPreamble),1,size(signal_rx,3));
%     transmitted_energy = vecnorm(signal,2,1).^2;%*vecnorm(H,2,1).^2;%*sum(vecnorm(H,2,1).^2,'all');           % energy of the signal
    noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
%     noise_var = noise_energy/(params.nActiveQ*(params.nData + params.nPreamble));   % variance of noise to be added
    noise_var = noise_energy/length(signal_tx1(:));
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
    
%     power = var(signal_tx1+signal_tx2);
%     noise_std = sqrt(power/(10^(SNR/10)));
    
    noise1 = noise_std.*(randn([size(signal_tx1),Nr])+1i*randn([size(signal_tx1),Nr]));      % noise
    noise2 = noise_var.*(rand([size(signal_tx1),Nr])+1i*rand([size(signal_tx1),Nr]))/2;      % noise
    
    % Noise added after the channel
    signal_rx = signal_rx1 + signal_rx2 +fft(noise1); %+ noise2;    % Add noise to signal
%     signal_rx = ifft(signal_rx1,[],1) + ifft(signal_rx2,[],1) +noise1;
%     signal_rx = fft(signal_rx,[],1);
    
end