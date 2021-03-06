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
   
function [Qam_symb] = modulation(params,bits,modulation)
    
%     % Number of bits per symbol
%     params.Nbps = params.modulation.params.Nbps;
Nsymb = size(bits,1)/params.Nbps; % Number of symbols
bits2 = reshape(bits,params.Nbps,Nsymb)';


switch modulation
    
    case 'qpsk'
        
        % Number of qam symbol
        N_qam_symbol = size(bits,1)/params.Nbps;
        % Generation of the qam modulation
        Qam_symb = zeros(N_qam_symbol,1);
        % Reshape the bits in (N_qam_symbol,params.Nbps)
        bits_tx = reshape(bits,params.Nbps,N_qam_symbol).';
    
        % In-Phase part
        params.NbpsI = params.Nbps/2;
        % Take the columns which correspond to the in-phase component
        bit_txI = bits_tx(:,1:params.NbpsI);
        % Gray to binary
        mapp_txI(:,1) = bit_txI(:,1);
        for ii = 2:params.NbpsI
            mapp_txI(:,ii) = xor(mapp_txI(:,ii-1) , bit_txI(:,ii)); 
        end

        % Binary to integer
        int_txI = bi2de(fliplr(mapp_txI));

        % Integer to symbol
        sigmaI = sqrt(sum(([0:2^params.NbpsI-1]-(2^params.NbpsI-1)/2).^2)/2^params.NbpsI); 
        symb_txI = 1/sigmaI/sqrt(2) * (int_txI - (2^params.NbpsI-1)/2);
        
        
        % Quadratic-Phase part
        params.NbpsQ = params.Nbps/2;
        bit_txQ = bits_tx(:,params.NbpsQ+1:end);
        
        % Gray to binary
        mapp_txQ(:,1) = bit_txQ(:,1);
        for ii = 2:params.NbpsQ
            mapp_txQ(:,ii) = xor( mapp_txQ(:,ii-1) , bit_txQ(:,ii) );
        end
        
        % Binary to integer
        int_txQ = bi2de(fliplr(mapp_txQ));
        
        % Integer to symbol
        sigmaQ = sqrt(sum(([0:2^params.NbpsQ-1]-(2^params.NbpsQ-1)/2).^2)/2^params.NbpsQ);
        symb_txQ = 1/sigmaQ/sqrt(2) * (int_txQ - (2^params.NbpsQ-1)/2);
        
        % Complex symbol
        Qam_symb = symb_txI + 1i*symb_txQ;
        
    case 'bpsk'
        % Gray to binary
        mapp_tx(:,1) = bits2(:,1);
        for ii = 2:params.Nbps
           mapp_tx(:,ii) = xor( mapp_tx(:,ii-1) , bits2(:,ii) ); 
        end

        % Binary to integer
        int_tx = bi2de(fliplr(mapp_tx)); %

        % Integer to symbol
        sigma = sqrt(sum(([0:2^params.Nbps-1]-(2^params.Nbps-1)/2).^2)/2^params.Nbps); 
        Qam_symb = 1/sigma * (int_tx - (2^params.Nbps-1)/2);
end
        
        
    % ---------------------------------------------------------------------
    % 'simple_modulation' just generates a certain amount of bits without 
    % any format and maps them into QAM symbols.
    % IMPORTANT!!: Comment the next line when trying your implementation
    % [bits,Qam_symb] = simple_modulation(params,Nbits);
end

