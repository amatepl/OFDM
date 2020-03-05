function [bits_pre,bits_mes,bits_pilot] = build_message(params,Nbits)

    %% -- Preamble
    % 1. construction of the preamble 2 OFDM symbols [preamble, preamble]
    %    without the inactive subcarrier.
    N_subcrr_act = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr;
    bits_pre= randi([0 1], N_subcrr_act * params.modulation.Nbps,1);
    % Two copies of the preamble follow eachother
    bits_pre = vertcat(bits_pre,bits_pre);
    
    % bits(1:N_subcrr_act * params.modulation.Nbps,1)=Preamble;   % Copy #1
    % bits(N_subcrr_act * params.modulation.Nbps+1:2*N_subcrr_act * params.modulation.Nbps,1)=Preamble;   % Copy #2
    
    %% -- Random message 
    % 2. Creation of the bit message without the inactive subcarrier and the pilots
    %    [message_1, message_2, ..., message_#ofdmsymbol]
    %N_subcrr_act = params.ofdm.N_subcrr - params.ofdm.N_inactive_subcrr - params.ofdm.N_pilots;
    bits_mes = randi([0 1], Nbits,1);
    
    %% -- Pilot
    % 3. Creation of the pilot [mes-pilot-a-pilot-ge]
    bits_pilot = randi([0 1],params.modulation.Nbps,1);
   
    
    