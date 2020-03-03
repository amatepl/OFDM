function [bits_pre,bits_mes,bits_pilot] = build_message(params,Nbits)

    % Preamble
    N_subcrr_act = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;
    bits_pre= randi([0 1], N_subcrr_act * params.modulation.Nbps,1);
    % Two copies of the preamble fallow eachother
    bits_pre = vertcat(bits_pre,bits_pre);
    
%     bits(1:N_subcrr_act * params.modulation.Nbps,1)=Preamble;   % Copy #1
%     bits(N_subcrr_act * params.modulation.Nbps+1:2*N_subcrr_act * params.modulation.Nbps,1)=Preamble;   % Copy #2
    
    % Randome message
    %N_subcrr_act = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr - params.ofdm.N_pilots;
    bits_mes = randi([0 1], Nbits,1);
    
    % Pilot
    bits_pilot = randi([0 1],params.modulation.Nbps,1);
    
    