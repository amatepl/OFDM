function STO = estimationSTO(params,signal_rx,maxiter)
    k=size(signal_rx,1); 
    STO=zeros(k,1);
    T = 1/params.B;
    N = params.Q+params.LCP;
    length_analyzed_signal = maxiter*(params.nPreamble+params.nData+params.N_zeros)*N; 
    
    for i=1:k
        signalrx=reshape(signal_rx(i,:),size(signal_rx,2),1);
        
        % Compute the correlation
        signalrx_mult = signalrx(1:length_analyzed_signal).*conj(signalrx(N+1:(length_analyzed_signal)+N));
        o = ones(N,1);
        An = conv(signalrx_mult,o,'valid');
        
        % Compute the norm of two following symbols
        norm = abs(signalrx(1:(length_analyzed_signal)+N)).^2;
        o = ones(2*N,1);
        norm = conv(norm,o,'valid');
        
        n = abs(An)./norm;
        
        STO(i) = find(n==max(n))-1;        
    end
end