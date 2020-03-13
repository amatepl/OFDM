function [STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx)
    k=size(signal_rx,1); 
    STO_estimated=zeros(k,1);
    CFO_estimated = zeros(k,1);
    T = 1/params.B;
    N = (params.Q+params.LCP);
    An = zeros(N,1);
    n = An;
    
%   Optimisation under work    
%     o = triu(tril(ones(size(signal_rx(N+1:(2*N)+N-1)))),-N+1);
%     o = o(:,1:N);
%     o = kron(o,ones(1,1,k));
% 
%     o = reshape(signal_rx.',[],1,k).*o;
%     
%     An = conj(signal_rx(1:2*N-1)*o);
    
    
    for i=1:k
        
        signalrx=reshape(signal_rx(i,:),1,size(signal_rx,2));
        % M is a matrix where every column is the same and we apply over it a
        % lower triangular transformation followed by a upper triangular one.
        
        
        
%         M = triu(tril(signalrx(N+1:(2*N)+N-1)'.*ones(size(signalrx(N+1:(2*N)+N-1)))),-N+1);
%         M = M(:,1:N);
%        
%         An = conj(signalrx(1:2*N-1)*M);
%     
%         M = triu(tril(signalrx(1:(2*N)+N-1)'.*ones(size(signalrx(1:(2*N)+N-1)))),-2*N +1);
%         M = M(:,1:N);
%     
%         n = abs(An)./(vecnorm(M).^2);
        
        for j=1:174080
            An(j) = signalrx(j:N+j-1)*(signalrx(N+j:(2*N)+j-1)');
            n(j) = norm(An(j))/norm(signalrx(j:(2*N)+j))^2;
        end
        figure;
        plot(n);
        
        STO_estimated(i) = find(n==max(n))-1;
        CFO_estimated(i) = angle(An(STO_estimated(i)+1))/(T*N);
        %CFO_estimated = 0;
        
    end
end