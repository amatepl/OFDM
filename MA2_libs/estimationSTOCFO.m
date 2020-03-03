function [STO_estimated, CFO_estimated] = estimationSTOCFO(params,signal_rx)
     
    T = 1/params.ofdm.B;
    % SFO & CFO Estimation & Correction
    

    N = params.ofdm.N_subcrr+params.ofdm.cp_L;
    %An = zeros(N,1);
    %n = An;
    
    % M is a matrix where every column is the same and we apply over it a
    % lower triangular transformation followed by a upper triangular one.
    M = triu(tril(signal_rx(N+1:(2*N)+N-1)'.*ones(size(signal_rx(N+1:(2*N)+N-1)))),-N+1);
    M = M(:,1:N);
    
    An = conj(signal_rx(1:2*N-1)*M);
    
    M = triu(tril(signal_rx(1:(2*N)+N-1)'.*ones(size(signal_rx(1:(2*N)+N-1)))),-2*N +1);
    M = M(:,1:N);
    
    n = abs(An)./(vecnorm(M).^2);
    
%     for i=1:N
%         An(i) = signal_rx(i:N+i-1)*(signal_rx(N+i:(2*N)+i-1)');
%         n(i) = norm(An(i))/norm(signal_rx(i:(2*N)+i-1))^2;
%     end
    STO_estimated = find(n==max(n))-1;
    CFO_estimated = -angle(An(STO_estimated+1))/(T*N);
    %CFO_estimated = 0;
end