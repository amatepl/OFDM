function [STO_estimated, CFO_estimated] = estimationSTOCFO(params,signal_rx)
    k=size(signal_rx,1); 
    STO_estimated=zeros(k,1);
    CFO_estimated = zeros(k,1);
    for i=1:k
        T = 1/params.ofdm.B;
        % SFO & CFO Estimation & Correction   
        
        N = params.ofdm.N_subcrr+params.ofdm.cp_L;
        An = zeros(N,1);
        n = An;
        
        % M is a matrix where every column is the same and we apply over it a
        % lower triangular transformation followed by a upper triangular one.
       % M = triu(tril(signal_rx(N+1:(2*N)+N-1,i)'.*ones(size(signal_rx(N+1:(2*N)+N-1,i)))),-N+1);
        %M = M(:,1:N);
        
        signalrx=reshape(signal_rx(i,:),1,size(signal_rx,2));
        %An = conj(signalrx(1:2*N-1)*M);
        
        
        %M = triu(tril(signal_rx(1:(2*N)+N-1,i)'.*ones(size(signal_rx(1:(2*N)+N-1,i)))),-2*N +1);
        %M = M(:,1:N);
        
        %n = abs(An)./(vecnorm(M).^2);
        
        
       for j=1:N
                An(j) = signalrx(1,j:N+j-1)*(signalrx(1,N+j:(2*N)+j-1)');
                n(j) = norm(An(j))/norm(signalrx(1,j:(2*N)+j-1))^2;
        end
        STO_estimated(i) = find(n==max(n))-1;
        CFO_estimated(i) = -angle(An(STO_estimated(i)+1))/(T*N);
        %CFO_estimated = 0;
    end
end