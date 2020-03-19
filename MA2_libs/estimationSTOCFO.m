function [STO_estimated, CFO_estimated] = estimationSTOCFO(params,signal_rx)
    k=size(signal_rx,1); 
    STO_estimated=zeros(k,1);
    CFO_estimated = zeros(k,1);
    
    T = 1/params.ofdm.B;
    N = params.ofdm.N_subcrr+params.ofdm.cp_L;
    
    length_analyzed_signal = 2*N*10; % Change this value to set the length of the signal you want to analyse 32*N is one frame.
    
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
        
        STO_estimated(i) = find(n==max(n))-1;
        CFO_estimated(i) = -angle(An(STO_estimated(i)+1))/(T*N);
        
        
%         signalrx=reshape(signal_rx(i,:),1,size(signal_rx,2));
%         % M is a matrix where every column is the same and we apply over it a
%         % lower triangular transformation followed by a upper triangular one.
%         M = triu(tril(signalrx(N+1:(2*N)+N-1)'.*ones(size(signalrx(N+1:(2*N)+N-1)))),-N+1);
%         M = M(:,1:N);
%        
%         An = conj(signalrx(1:2*N-1)*M);
%     
%         M = triu(tril(signalrx(1:(2*N)+N-1)'.*ones(size(signalrx(1:(2*N)+N-1)))),-2*N +1);
%         M = M(:,1:N);
%     
%         n = abs(An)./(vecnorm(M).^2);
%         
%         STO_estimated(i) = find(n==max(n))-1;
%         CFO_estimated(i) = -angle(An(STO_estimated(i)+1))/(T*N);
%         %CFO_estimated = 0;
    end
end