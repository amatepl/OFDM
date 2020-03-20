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
   
    H = permute(H,[2 3 1]); % Tranforms H into a 2048x1x4 matrix
    
    signal_rx = signal_tx.*H;   % Multiplication in frequency domain of the channel and the signal
    
     % 2. Noise computation
    transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
    noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
    noise_var = noise_energy/(length(signal_tx(:)));   % variance of noise to be added
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
    noise = noise_std*(randn(size(signal_tx))+1i*randn(size(signal_tx)));      % noise
    
    signal_rx = signal_rx +noise;    % Add noise to signal
    
end