function [STO_estimated, CFO_estimated] = estimationSTOCFO(params,signal_rx)
    N = params.ofdm.N_subcrr+params.ofdm.cp_L;
    An = zeros(N,1);
    n = An;
    for i=1:N
        An(i) = signal_rx(i:N+i-1)*(signal_rx(N+i:(2*N)+i-1)');
        n(i) = norm(An(i))/norm(signal_rx(i:(2*N)+i-1))^2;
    end
    STO_estimated = find(n==max(n))-1;
    CFO_estimated = 0;
end