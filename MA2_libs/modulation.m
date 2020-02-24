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
%           Nbits   : Number of bits. dim=(1,1)
% Outputs:
%           bits    : randomly generated bits. dim=(Nbits,1)
%           Qam_symb: result of mapping the bits into the 
%                     QAM constellation. dim=(N_qam_symb,1)
%
   
function [bits,Qam_symb] = modulation(params,Nbits)
    % **** YOUR CODE GOES HERE!!
    
    
    
    
    % ---------------------------------------------------------------------
    % 'simple_modulation' just generates a certain amount of bits without 
    % any format and maps them into QAM symbols.
    % IMPORTANT!!: Comment the next line when trying your implementation
    [bits,Qam_symb] = simple_modulation(params,Nbits);
end