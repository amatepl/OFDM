function [bits_pre,bits_mes,bits_pilot] = build_message_test(params,Nbits)

    % Preamble
    N_subcrr_act = params.nActiveQ;
    bits_pre= randi([0 1], N_subcrr_act * params.Nbps,1);
    % Two copies of the preamble fallow eachother
    bits_pre = vertcat(bits_pre,bits_pre);
    
    % Random message
    bits_mes = randi([0 1], Nbits,1);
    
    % Pilot
    bits_pilot = randi([0 1],params.Nbps,1);
    