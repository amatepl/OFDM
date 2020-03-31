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
    a = ones(5,Nr); 
    impulse_response(1,:) = a(1,:);
    % MPC
%     a(1,:) = a(1,:);% attenuation factor
%     a(2,:) = a(2,:).*rand(1,Nr)*0.8;                               % Constant attenuation
%     a(3,:) = a(3,:).*rand(1,Nr)*0.7;
%     a(4,:) = a(4,:).*rand(1,Nr)*0.5; 
%     a(5,:) = a(5,:).*rand(1,Nr)*0.3;
%     phi = 2*pi*randn(4,Nr);                             % random phase
%     phi = mod(phi,2*pi);                               % map it to the range [0,2*pi]
%     for i = 1:Nr
%         impulse_response(2+randi(20),i) = a(2,i).*exp(1i*phi(1,i));
%         impulse_response(21+randi(20),i) = a(3,i).*exp(1i*phi(2,i));
%         impulse_response(41+randi(20),i) = a(4,i).*exp(1i*phi(3,i));
%         impulse_response(61+randi(20),i) = a(5,i).*exp(1i*phi(4,i));
%     end
    %impulse_response = [zeros(STO,1); impulse_response(1:end-STO)];
    h = impulse_response(1:params.Q);
    H1 = fft(h,params.Q);
    H = zeros(params.Q,1);
    H(params.ActiveQIndex) = H1(params.ActiveQIndex);
    
%     figure;
%     stem(0:params.Q-1,abs(h));
%     xlim([0 20])
%     xlabel("time");
%     ylabel("|h(t)|");
%     grid on;
%     title("Channel impulse response in time domain with STO delay");
    
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
    Nbits = (params.nPreamble+params.nData)*(params.nActiveQ)*params.Nbps;
    signalEnergy = norm(signal_tx)^2;
    Energybit = signalEnergy/Nbits;
    EbNoLin = 2*10^(SNR/10);
    No = Energybit/EbNoLin;
    noise = sqrt(No/2)*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));

    % 4. STO addition
    signal_rx = [signal_rx(:,1:STO), signal_rx(:,1:end-STO)]; 
    
    % 5. Frequency offset (CFO)
    % delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
    delta_w = params.Fc*CFO;
    T = 1/params.B;
    n = 1:1:size(signal_rx,2);
    phi = exp(1i*delta_w*T*n);
    signal_rx = signal_rx.*phi;
    
    signal_rx = signal_rx+noise.'; 
    
    % 6. Matched filter + MMSE equalizer
    % impulse_response = [zeros(STO,1); impulse_response(1:end-STO)];
    % impulse_matrix = convolutionMatrix(impulse_response);
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