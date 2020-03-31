
function Qsymb_rx = simpleReceiver(params,signal_rx,preamble)
    signal_rx = reshape(signal_rx.',[],params.nData+params.nPreamble);
    signal_rx=signal_rx(params.LCP+1:end,:);
    signal_rx=fft(signal_rx(:,1:end),params.Q);
    signal_rx = signal_rx(params.ActiveQIndex,:);
    
    a=conj(preamble).*signal_rx(:,2);
    ht = ifft(a,params.nActiveQ,1);
    Ht=fft(ht,params.nActiveQ,1); 
    
%     signal_rx=signal_rx.*conj(Ht);
    
    signal_rx = signal_rx(:,params.nPreamble+1:end);
    
    symb_rx1 = signal_rx(1:size(signal_rx,1)/2,:);
    symb_rx2 = signal_rx(size(signal_rx,1)/2+1:end,:);
    if params.N_pilots > 0
        symb_rx1 = reshape(symb_rx1,[],params.nData,params.N_pilots/2);
        symb_rx2 = reshape(symb_rx2,[],params.nData,params.N_pilots/2);
        symb_rx1 = symb_rx1(2:end,:,:);
        symb_rx2 = symb_rx2(1:end-1,:,:);
        symb_rx1 = reshape(symb_rx1,[],params.nData);
        symb_rx2 = reshape(symb_rx2,[],params.nData);
    end
    Qsymb_rx = [symb_rx1;symb_rx2];
    Qsymb_rx = reshape(Qsymb_rx,[],1);
end

