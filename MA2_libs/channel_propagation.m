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

function signal_rx = channel_propagation(params,signal_tx,SNR)
    % **** YOUR CODE GOES HERE!!
    % Multipath component:
    Lcp = params.ofdm.cp_L;
    Q = params.ofdm.N_subcrr;
    signal_tx_col = reshape(signal_tx,Lcp+Q,[]);
    impulse_response = zeros(Lcp+Q,1);
    impulse_response(1) = 1;
    a = (rand(1,Lcp+Q-1)).';
    phi = 2*pi*(randn(1,Lcp+Q-1)).';
    % map it to the range [0,2*pi]
    phi = 0;%mod(phi,2*pi);
    impulse_response(2:end) = a.*exp(1i*phi);
    impulse_matrix = convolutionMatrix(impulse_response);
    signal_rx = impulse_matrix*signal_tx_col;
    signal_rx = reshape(signal_rx.',[],size(signal_rx,1)*size(signal_rx,2));
    
    transmitted_energy = norm(signal_tx(:))^2;           % energy of the signal
    noise_energy = transmitted_energy/(10^(SNR/10));     % energy of noise
    noise_var = noise_energy/(length(signal_tx(:))-1);   % variance of noise to be added
    noise_std = sqrt(noise_var/2);                       % std. deviation of noise to be added
    noise = noise_std*(randn(length(signal_tx),1)+1i*randn(length(signal_tx),1));      % noise
    signal_rx = (signal_rx).'+noise; 
%     % Matched filter + MMSE equalizer
%     Nb = 2 * params.ofdm.N_subcrr * params.modulation.Nbps;
%     No = noise_energy/(Nb*2);
%     signal_rx = ((2*No/var(signal_tx)+impulse_matrix'*impulse_matrix)\impulse_matrix'*signal_rx).';
   
    
    % ---------------------------------------------------------------------
    % 'simple_channel_AWGN': Implements a simple AWGN channel. No
    % STO, Multipath, etc.
    % IMPORTANT!!: Comment the next line when trying your implementation    
    % signal_rx = simple_channel_AWGN(params,signal_tx,SNR);
end