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
    % **** YOUR CODE GOES HERE!!
%     % 1. Some constants:
%     lambda = physconst('LightSpeed')/params.Fc;
%     beta = 2*pi/lambda;
%     Lcp = params.LCP;
%     Q = params.Q;
%     
%     % 2. Multipath component
%     % Reshape the signal tx in a matrix 
%     signal_tx_col = reshape(signal_tx,Lcp+Q,[]);
%     % Creation of the impulse response matrix
%     impulse_response = zeros(Lcp+Q,Nr);
%     a = ones(2,Nr);                                   % attenuation factor
%     a(2,:) = 0.5*a(2,:);
%     phi = 2*pi*rand(2,Nr);                            % random phase
%     phi = mod(phi,2*pi);                              % map it to the range [0,2*pi]
%     z = lambda*(linspace(0,Nr-1,Nr)/2-(Nr-1)/4);      % z coordinate in the local zone
%     u = ones(2,Nr);                                   % angular transmission
%     u(1,:) = beta*cos(pi/6)*u(1,:);
%     u(2,:) = beta*cos(2*pi/3)*u(2,:);
%     impulse_response(1,:) = a(1,:);
%     %impulse_response(1,:) = a(1,:).*exp(1i*phi(1,:)).*exp(1i*u(1,:).*z);
%     impulse_response(4,:) = a(2,:).*exp(1i*phi(2,:)).*exp(1i*u(2,:).*z);
%     impulse_matrix = convolutionMatrix(impulse_response,Nr);
%     
%     % Convolution of the signal with the impulse matrix
%     signal_rx_SIMO = zeros(size(signal_tx_col,1),size(signal_tx_col,2));
%     signal_rx = zeros(Nr,size(signal_rx_SIMO,1)*size(signal_rx_SIMO,2));
%     
%     for i = 1:Nr
%         signal_rx_SIMO = impulse_matrix(:,:,i)*signal_tx_col;
%         signal_rx(i,:) = reshape(signal_rx_SIMO,size(signal_rx_SIMO,1)*size(signal_rx_SIMO,2),1).';
%     end
%     
%     clear signal_rx_SIMO;
%     clear impulse_response;


%     figure, hold on;
%     plot(abs(H(1,:).'));

    h = ifft(H,[],2);
    signal_rx = zeros(Nr,length(signal_tx));
%     signal_rx = conv2(1,signal_tx,h);
    
    for i = 1:Nr
        signal_rx(i,:) = conv(signal_tx,h(i,:),'same');
    end


%     signal_rx = signal_tx.*ones(Nr,length(signal_tx));
    
    % 3. Noise computation
    transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
    noise_energy = transmitted_energy/(2*10^(SNR/10));     % energy of noise
    noise_var = noise_energy/(length(signal_tx(:)));   % variance of noise to be added
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
    noise = noise_std*(randn(length(signal_tx),Nr)+1i*randn(length(signal_tx),Nr));      % noise
    
    
    signal_rx = signal_rx +noise.'; 
    
    % Matched filter + MMSE equalizer
%     impulse_response = [zeros(STO,1); impulse_response(1:end-STO)];
%     impulse_matrix = convolutionMatrix(impulse_response);
%     Nb = 2 * params.ofdm.N_subcrr * params.modulation.Nbps;
%     No = noise_energy/(Nb*2);
%     signal_rx_col = reshape(signal_rx,Lcp+Q,[]);
%     var1 = var(signal_rx);
%     signal_rx_col = (2*No/var1+impulse_matrix'*impulse_matrix)\impulse_matrix'*signal_rx_col;
%     signal_rx = reshape(signal_rx_col,size(signal_rx,1)*size(signal_rx,2),1).';
   
    
    % ---------------------------------------------------------------------
    % 'simple_channel_AWGN': Implements a simple AWGN channel. No
    % STO, Multipath, etc.
    % IMPORTANT!!: Comment the next line when trying your implementation    
    % signal_rx = simple_channel_AWGN(params,signal_tx,SNR);
end