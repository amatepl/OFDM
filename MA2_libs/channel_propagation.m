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
    % signal_rx = awgn(signal_tx,SNR,'measured');
    sigEner = norm(signal_tx(:))^2;                    % energy of the signal
    noiseEner = sigEner/(10^(SNR/10));        % energy of noise to be added
    noiseVar = noiseEner/(length(signal_tx(:)-1));     % variance of noise to be added
    noiseStd = sqrt(noiseVar);                   % std. deviation of noise to be added
    noise = noiseStd/2*(randn(size(signal_tx))+1i*randn(size(signal_tx)));           % noise
    signal_rx = signal_tx+noise; 
    
    
    
    % ---------------------------------------------------------------------
    % 'simple_channel_AWGN': Implements a simple AWGN channel. No
    % STO, Multipath, etc.
    % IMPORTANT!!: Comment the next line when trying your implementation    
    % signal_rx = simple_channel_AWGN(params,signal_tx,SNR);
end