% Opera - ULB
% February 2020
%
% ToDo: 
%       - Define the format of your frame: preamble and data
%       - Generate the bits for your frame.
%       - Map the generated bits into QAM symbols.
%
% Inputs: 
%           params  : MA2 parameters. 
%           symb_rx : received qam symbols. dim=(N_qam_symb,1)
% Outputs:
%           bit_rx  : received bits. dim=(Nbits,1)
%
   
function bits_rx = demodulation(params,symb_rx)
    % **** YOUR CODE GOES HERE!!
    % *** If you feel like don't touch this part but ensure the input data
    % *** is in the appropriate format ;)
    
    
    
    
    % ---------------------------------------------------------------------
    % 'demapping' Simple mapping of qam symbols to bits
    % IMPORTANT!!: Comment the next line when trying your implementation
    bits_rx = demapping(symb_rx,params.modulation.Nbps);
end