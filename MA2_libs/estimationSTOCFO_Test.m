function [STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx)
    k=size(signal_rx,1); 
    STO_estimated=zeros(k,1);
    CFO_estimated = zeros(k,1);
    for i=1:k
        T = 1/params.B;
        N = params.Q+params.LCP;
        signalrx=reshape(signal_rx(i,:),1,size(signal_rx,2));
        % M is a matrix where every column is the same and we apply over it a
        % lower triangular transformation followed by a upper triangular one.
        M = triu(tril(signalrx(N+1:(2*N)+N-1)'.*ones(size(signalrx(N+1:(2*N)+N-1)))),-N+1);
        M = M(:,1:N);
       
        An = conj(signalrx(1:2*N-1)*M);
    
        M = triu(tril(signalrx(1:(2*N)+N-1)'.*ones(size(signalrx(1:(2*N)+N-1)))),-2*N +1);
        M = M(:,1:N);
    
        n = abs(An)./(vecnorm(M).^2);
        
        STO_estimated(i) = find(n==max(n))-1;
        CFO_estimated(i) = -angle(An(STO_estimated(i)+1))/(T*N);
        %CFO_estimated = 0;
    end
end