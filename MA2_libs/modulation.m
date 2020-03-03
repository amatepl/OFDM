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
   
function [Qam_symb] = modulation(Nbps,bits)
    
%     % Number of bits per symbol
%     Nbps = params.modulation.Nbps;

    % Number of qam symbol
    N_qam_symbol = size(bits,1)/Nbps;
    % Generation of the qam modulation
    Qam_symb = zeros(N_qam_symbol,1);
    % Reshape the bits in (N_qam_symbol,Nbps)
    bits_tx = reshape(bits,Nbps,N_qam_symbol).';
    
    % In-Phase part
    NbpsI = Nbps/2;
    % Take the columns which correspond to the in-phase component
    bit_txI = bits_tx(:,1:NbpsI);
    % Gray to binary
    mapp_txI(:,1) = bit_txI(:,1);
    for ii = 2:NbpsI
       mapp_txI(:,ii) = xor(mapp_txI(:,ii-1) , bit_txI(:,ii)); 
    end

    % Binary to integer
    int_txI = bi2de(fliplr(mapp_txI));

    % Integer to symbol
    sigmaI = sqrt(sum(([0:2^NbpsI-1]-(2^NbpsI-1)/2).^2)/2^NbpsI); 
    symb_txI = 1/sigmaI/sqrt(2) * (int_txI - (2^NbpsI-1)/2);
        
        
    % Quadratic-Phase part
    NbpsQ = Nbps/2;
    bit_txQ = bits_tx(:,NbpsQ+1:end);
        
    % Gray to binary
    mapp_txQ(:,1) = bit_txQ(:,1);
    for ii = 2:NbpsQ
       mapp_txQ(:,ii) = xor( mapp_txQ(:,ii-1) , bit_txQ(:,ii) ); 
    end

    % Binary to integer
    int_txQ = bi2de(fliplr(mapp_txQ));

    % Integer to symbol
    sigmaQ = sqrt(sum(([0:2^NbpsQ-1]-(2^NbpsQ-1)/2).^2)/2^NbpsQ); 
    symb_txQ = 1/sigmaQ/sqrt(2) * (int_txQ - (2^NbpsQ-1)/2);
       
    % Complex symbol
    Qam_symb = symb_txI + 1i*symb_txQ;
    % ---------------------------------------------------------------------
    % 'simple_modulation' just generates a certain amount of bits without 
    % any format and maps them into QAM symbols.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % [bits,Qam_symb] = simple_modulation(params,Nbits);
end

