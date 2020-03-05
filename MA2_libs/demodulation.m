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
   
function bits_rx = demodulation(params,symb_rx, modulation)
   
% INPUTS:
% - symb_rx : vector of input symbols (variance 1)
% - Nbps : number of bits per symbol
% - modulation : 'pam' or 'qam'
%
% OUTPUTS:
% - bit_rx : vector of ouput bits

    switch modulation
        
        case 'qpsk'
            
            Nsymb = size(symb_rx,1); % Number of symbols
            % Number of bits per symbol
%             Nbps = params.modulation.Nbps;
            Nbps = 2;
            % REAL PART
            NbpsI = Nbps/2;
            symb_rxI = real(symb_rx);
            
            % Symbol to integer
            sigmaI = sqrt(sum(([0:2^NbpsI-1]-(2^NbpsI-1)/2).^2)/2^NbpsI);
            int_rxI = sigmaI * sqrt(2) * symb_rxI + (2^NbpsI-1)/2;
            
            % Integer detection
            int_detI = round(int_rxI);
            int_detI(int_detI<0) = 0;
            int_detI(int_detI>2^NbpsI-1) = 2^NbpsI-1;
            
            % Integer to binary
            mapp_rxI  = fliplr(de2bi(int_detI));
            
            % Binary to gray
            bit_rxI(:,1) = mapp_rxI(:,1);
            for ii = 2:NbpsI
                bit_rxI(:,ii) = xor( mapp_rxI(:,ii-1) , mapp_rxI(:,ii) );
            end
            
            
            % IMAGINARY PART
            NbpsQ = Nbps/2;
            symb_rxQ = imag(symb_rx);
            
            % Symbol to integer
            sigmaQ = sqrt(sum(([0:2^NbpsQ-1]-(2^NbpsQ-1)/2).^2)/2^NbpsQ);
            int_rxQ = sigmaQ * sqrt(2) * symb_rxQ + (2^NbpsQ-1)/2;
            
            % Integer detection
            int_detQ = round(int_rxQ);
            int_detQ(int_detQ<0) = 0;
            int_detQ(int_detQ>2^NbpsI-1) = 2^NbpsQ-1;
            
            % Integer to binary
            mapp_rxQ  = fliplr(de2bi(int_detQ));
            
            % Binary to gray
            bit_rxQ(:,1) = mapp_rxQ(:,1);
            for ii = 2:NbpsQ
                bit_rxQ(:,ii) = xor( mapp_rxQ(:,ii-1) , mapp_rxQ(:,ii) );
            end
            
            % BIT CONCATENATION
            bits_rx = reshape([bit_rxI,bit_rxQ]',Nsymb*Nbps,1);
            
        case 'bpsk'
            Nbps = 1;
            % Symbol to integer
            sigma = sqrt(sum(([0:2^Nbps-1]-(2^Nbps-1)/2).^2)/2^Nbps); 
            int_rx = sigma * symb_rx + (2^Nbps-1)/2;
            
            % Integer detection
            
            int_det = round(real(int_rx));
            %         int_det(find(int_det<0)) = 0;
            %         int_det(find(int_det>2^Nbps-1)) = 2^Nbps-1;
            int_det(int_det<0) = 0;
            int_det(int_det>2^Nbps-1) = 2^Nbps-1;
            
            % Integer to binary
            mapp_rx  = fliplr(de2bi(int_det));
            
            % Binary to gray
            bit_rx2(:,1) = mapp_rx(:,1);
            for ii = 2:Nbps
                bit_rx2(:,ii) = xor( mapp_rx(:,ii-1) , mapp_rx(:,ii) );
            end
            Nsymb = size(symb_rx,1); % Number of symbols
            bits_rx = reshape(bit_rx2',Nsymb*Nbps,1);
            
    end
            
            
    
    % ---------------------------------------------------------------------
    % 'demapping' Simple mapping of qam symbols to bits
    % IMPORTANT!!: Comment the next line when trying your implementation
    %bits_rx = demapping(symb_rx,params.modulation.Nbps);
end