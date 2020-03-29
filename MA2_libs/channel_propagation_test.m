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

function signal_rx = channel_propagation_test(params,signal_tx,SNR,STO,CFO,Nr)
    % **** YOUR CODE GOES HERE!!
    
    % 1. Some constants:
    lambda = physconst('LightSpeed')/params.Fc;
    beta = 2*pi/lambda;
    Lcp = params.LCP;
    Q = params.Q;
    
    % 2. Multipath component
    % Reshape the signal tx in a matrix 
    signal_tx_col = reshape(signal_tx,Lcp+Q,[]);
    % Creation of the impulse response matrix
    impulse_response = zeros(Lcp+Q,Nr);
    a = ones(2,Nr);                                    % attenuation factor
    a(2,:) = 0.4*a(2,:);                               % Constant attenuation
    phi = 2*pi*randn(1,Nr);                             % random phase
    phi = mod(phi,2*pi);                               % map it to the range [0,2*pi]
    impulse_response(1,:) = a(1,:);
    impulse_response(2,:) = a(2,:).*exp(1i*phi(1,:));
    h = impulse_response(1:params.Q);
    H1 = fft(h,params.Q);
    H = zeros(params.Q,1);
    H(params.ActiveQIndex) = H1(params.ActiveQIndex);
    
%     figure;
%     subplot(2,1,1);
%     stem(0:params.Q-1,abs(h));
%     xlim([0 10])
%     xlabel("time");
%     ylabel("|h(t)|");
%     grid on;
%     subplot(2,1,2);
%     plot(abs(fftshift(H)));
%     xlabel("frequency (bins)");
%     ylabel("|H(f)|");
%     grid on;
%     sgtitle("Channel impulse response in time and frequency domain");
    
    
    impulse_matrix = convolutionMatrix(impulse_response,Nr);
    
    % Convolution of the signal with the impulse matrix
    signal_rx_SIMO = zeros(size(signal_tx_col,1),size(signal_tx_col,2));
    signal_rx = zeros(Nr,size(signal_rx_SIMO,1)*size(signal_rx_SIMO,2));
    
    for i = 1:Nr
        signal_rx_SIMO = impulse_matrix(:,:,i)*signal_tx_col;
        signal_rx(i,:) = reshape(signal_rx_SIMO,size(signal_rx_SIMO,1)*size(signal_rx_SIMO,2),1).';
    end
    
    clear signal_rx_SIMO;
    clear impulse_response;
    
    % 3. Noise addition

%     transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
%     noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
%     noise_var = noise_energy/(length(signal_tx(:)));   % variance of noise to be added
%     noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
%     noise = noise_std*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));      % noise
    

%     transmitted_energy = trapz(abs(signal_tx)).^2;           % energy of the signal
%     noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
%     noise_var = noise_energy/(length(signal_tx));   % variance of noise to be added
%     noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
%     noise = noise_std*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));      % noise
    Nbits = (params.nPreamble+params.nData)*(params.nActiveQ)*params.Nbps;
    signalEnergy = norm(signal_tx)^2;
    Energybit = signalEnergy/Nbits;
    EbNoLin = 10^(SNR/10);
    No = Energybit/EbNoLin;
    noise = sqrt(No/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));

    % 4. STO addition
%     signal_rx = [zeros(Nr,STO), signal_rx(:,1:end-STO)]; 
    
    % 5. Frequency offset (CFO)
    % delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
%     delta_w = params.Fc*CFO*(2*pi);
%     T = 1/params.B;
%     n = 1:1:size(signal_rx,2);
%     phi = exp(1i*delta_w*T*n);
%     signal_rx = signal_rx.*phi;
    
    signal_rx = signal_rx;%+noise.'; 
    
    % 6. Matched filter + MMSE equalizer
    % impulse_response = [zeros(STO,1); impulse_response(1:end-STO)];
    % impulse_matrix = convolutionMatrix(impulse_response);
    % Nb = 2 * params.ofdm.N_subcrr * params.modulation.Nbps;
    % No = noise_energy/(Nb*2);
    % signal_rx_col = reshape(signal_rx,Lcp+Q,[]);
    % var1 = var(signal_rx);
    % signal_rx_col = (2*No/var1+impulse_matrix'*impulse_matrix)\impulse_matrix'*signal_rx_col;
    % signal_rx = reshape(signal_rx_col,size(signal_rx,1)*size(signal_rx,2),1).';
   
    
    % ---------------------------------------------------------------------
    % 'simple_channel_AWGN': Implements a simple AWGN channel. No
    % STO, Multipath, etc.
    % IMPORTANT!!: Comment the next line when trying your implementation    
    % signal_rx = simple_channel_AWGN(params,signal_tx,SNR);
end