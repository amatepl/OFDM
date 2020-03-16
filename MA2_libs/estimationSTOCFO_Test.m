function [STO_estimated, CFO_estimated] = estimationSTOCFO_Test(params,signal_rx)
    k=size(signal_rx,1); 
    STO_estimated=zeros(k,1);
    CFO_estimated = zeros(k,1);
    T = 1/params.B;
    N = params.Q+params.LCP;
    length_analyzed_signal = 2*32*N;
    
%      N = (params.nData + params.nPreamble )*(params.Q+params.LCP);
    
%   Optimisation under work    
%     o = triu(tril(ones(size(signal_rx(N+1:(2*N)+N-1)))),-N+1);
%     o = o(:,1:N);
%     o = kron(o,ones(1,1,k));
% 
%     o = reshape(signal_rx.',[],1,k).*o;
%     
%     An = conj(signal_rx(1:2*N-1)*o);

    m = 2;
    
    
    
    for i=1:k
        
        signalrx=reshape(signal_rx(i,:),size(signal_rx,2),1);
        
        
%         signalrx = signalrx(1:sig_to_analyse).*signalrx(1+N:sig_to_analyse+N);
%         signalrx1 = signalrx(1:m*N-1).*conj(signalrx(N+1:(m*N)+N-1));
        
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
        
%         % M is a matrix where every column is the same and we apply over it a
%         % lower triangular transformation followed by a upper triangular one.
%             
%         
% %         M = triu(tril(ones(size(signalrx(N+1:(m*N)+N-1),1))),-N+1);
%         
% %         signalrx =  signalrx(N+1:(m*N)+N-1)';
%         
% %         M = signalrx(N+1:(m*N)+N-1).*M;
%         
%         M = triu(tril(signalrx(N+1:(m*N)+N-1).*ones(size(signalrx(N+1:(m*N)+N-1),1))),-N+1);
%         M = M(:,1:N);
%        
% %         An = conj(signalrx(1:m*N-1)*M);
% 
%         An = (signalrx(1:m*N-1).')*M;
%     
%         M = triu(tril(signalrx(1:(m*N)+N-1).*ones(size(signalrx(1:(m*N)+N-1),1))),-2*N +1);
%         M = M(:,1:N);
%     
%         n = abs(An)./(vecnorm(M).^2);
%         
%         STO_estimated(i) = find(n==max(n))-1;
%         CFO_estimated(i) = -angle(An(STO_estimated(i)+1))/(T*N);
%         %CFO_estimated = 0;
        
    end
end