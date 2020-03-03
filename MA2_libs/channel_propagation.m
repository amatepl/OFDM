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

function signal_rx = channel_propagation(params,signal_tx,SNR,STO,CFO)
    % **** YOUR CODE GOES HERE!!
    % Multipath component:
    Lcp = params.ofdm.cp_L;
    Q = params.ofdm.N_subcrr;
    signal_tx_col = reshape(signal_tx,Lcp+Q,[]);
    impulse_response = zeros(Lcp+Q,1);
    impulse_response(1) = 1;
    %a = (rand(1,Lcp+Q-1)).';
    a = rand(1,1);
    %phi = 2*pi*(randn(1,Lcp+Q-1)).';
    phi = 2*pi*randn(1,1);
    % phi = 0;
    % map it to the range [0,2*pi]
    phi = mod(phi,2*pi);
    %impulse_response(randi([1,Lcp+Q])) = a*exp(1i*phi);
    impulse_response(2) = a*exp(1i*phi);
    impulse_matrix = convolutionMatrix(impulse_response);
    signal_rx = impulse_matrix*signal_tx_col;
    signal_rx = reshape(signal_rx,size(signal_rx,1)*size(signal_rx,2),1).';
    
    % Noise
    transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
    noise_energy = transmitted_energy/(10^(SNR/10));     % energy of noise
    noise_var = noise_energy/(length(signal_tx(:))-1);   % variance of noise to be added
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
    noise = noise_std*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));      % noise
    
    % STO
    signal_rx = signal_tx;
    signal_rx = [zeros(1,STO), signal_rx(1:end-STO)];
    
    % Frequency offset (CFO)
    
    % delta_w is usually in the range [-40ppm, 40ppm] Source: Wikipedia
    delta_w = params.ofdm.f_dc*CFO/(2*pi);
    T = 1/params.ofdm.B;
    n = 1:1:size(signal_rx,2);
    phi = exp(1i*delta_w*T*n);
    signal_rx = signal_rx.*phi;
    
    signal_rx = signal_rx+noise.'; 
    
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