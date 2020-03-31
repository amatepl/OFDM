clear, close all, clc;

cfg = load('MA2_lab_parameters.mat');   % load configFile
SNR_list = cfg.params.SNR_list;                    % get the set of parameters

addpath(genpath('Results'));           % add data

%% Paramters to set ------------------------------
Htype = 'LOS';          % Can be LOS or NLOS
System = 'MIMO';        % SIMO or MIMO
%% -----------------------------------------------
Nreceivers = 4;
Nbps = 2;

disp('$$ Displaying results:');
figure;
% ber theoretical
ber_theo = berawgn(SNR_list,'qam',2^(Nbps));
semilogy(SNR_list,ber_theo,'--');
grid on;hold on;

switch System
    case 'SIMO'
        title(join([Htype,' ',System,' system']));
%         for Nr = 1:Nreceivers
%             BER_i = load(join(['SIMO_',Htype,'_',num2str(Nr),'.mat']));
%             BER_i = BER_i.BER_i;
%             semilogy(SNR_list,mean(BER_i,1));
%         end
        
        BER_i = load(join(['SIMO_',Htype,'_','1','.mat']));
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1),'-o');
        
        BER_i = load(join(['SIMO_',Htype,'_','2','.mat']));
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1),'-+');
        
        BER_i = load(join(['SIMO_',Htype,'_','3','.mat']));
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1),'-*');
        
        BER_i = load(join(['SIMO_',Htype,'_','4','.mat']));
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1),'-^','DisplayName','M_T = 1, M_R = 4');
        
        legend('theoretical','M_T = 1, M_R = 1','M_T = 1, M_R = 2','M_T = 1, M_R = 3','M_T = 1, M_R = 4');
    case 'MIMO'
        title(join([Htype,' NLOS and LOS systems']));
        
        BER_i = load('MIMO_LOS_4.mat');
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1));
        
        title(join([Htype,'NLOS and LOS systems']));
        BER_i = load('MIMO_NLOS_4.mat');
        BER_i = BER_i.BER_i;
        semilogy(SNR_list,mean(BER_i,1));
        
        legend('theoretical','LOS','NLOS');
end


% legend('myBER','theoretical');
xlabel('SNR dB');ylabel('Probability of error');
xlim([-5 20]);
ylim([10^(-6) 1]);


