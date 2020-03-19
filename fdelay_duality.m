 clc; close all; clear;
 
 params = load('TestParam.mat').TestParam;

 H_tensor=load('H_tensor.mat').H_tensor;
 H_tensor_nlos=load('H_tensor_nlos.mat').H_tensor_nlos;
 
 %% fit rice distrib to channel frequency response for chosen antenna and chosen frame
 x_val=0:0.01:5; %values of |h| where the pdf is evaluated
 ylos=zeros(size(x_val,2),4); ynlos=zeros(size(x_val,2),4);
 Klos=zeros(4); Knlos=zeros(4);
 for i=1:4
    p1 = fitdist(abs(H_tensor(:,10,i)),'Rician');
    p2 = fitdist(abs(H_tensor_nlos(:,10,i)),'Rician');
    ylos(:,i)=pdf(p1,x_val);
    ynlos(:,i)=pdf(p2,x_val);
    Klos(i)=10*log10(p1.s^2/2*p1.sigma); %rice factor
    Knlos(i)=10*log10(p2.s^2/2*p1.sigma);
 end
 %histograms
 figure;
 for i=1:4
     subplot(2,2,i);
     plot(x_val,ylos(:,i));
     hold on;
     histogram(abs(H_tensor(:,10,i)),'BinMethod','auto','Normalization','pdf','DisplayStyle','stairs','LineStyle','-.');
     hold on;
     grid on;
     legend('fit', 'histogram');
     ylabel('PDF');
     xlabel('|h|');
     title(sprintf('antenna %0.5g',i));    
 end   
 sgtitle('Rice fitting channel frequency response LOS');
 
 figure;
  for i=1:4
     subplot(2,2,i);
     plot(x_val,ynlos(:,i));
     hold on;
     histogram(abs(H_tensor_nlos(:,10,i)),'BinMethod','auto','Normalization','pdf','DisplayStyle','stairs','LineStyle','-.');
     hold on;
     grid on;
     legend('fit', 'histogram');
     ylabel('PDF');
     xlabel('|h|');
     title(sprintf('antenna %0.5g',i));    
 end   
 sgtitle('Rice fitting channel frequency response NLOS');
 % Rice fits
 Legend=cell(1,4);
 figure;
 for i=1:4
    plot(x_val,ylos(:,i));
    hold on;
    Legend{i}=strcat(sprintf('antenna %0.5g: KdB =',i), num2str(Klos(i)));
 end
 grid on;
 xlabel('|h|');
 ylabel('PDF');
 legend(Legend);
 title('Rice distribution fitting on channel frequency response LOS');
 
 Legend=cell(1,4);
 figure;
 for i=1:4
    plot(x_val,ynlos(:,i));
    hold on;
    Legend{i}=strcat(sprintf('antenna %0.5g: KdB =',i), num2str(Knlos(i)));
 end
 grid on;
 xlabel('|h|');
 ylabel('PDF');
 legend(Legend);
 title('Rice distribution fitting on channel frequency response NLOS'); 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %%  Rice distrib fitting to channel impulse response taps LOS
 % taps distributions of CIR's
 
 ht = ifft(H_tensor(:,:,4));
 %ht=ht(1:256,:);
 
 htt=(abs(ht)).^2;
 avg_power_taps = 10*log10(mean(htt,2));
 strongest_tap = max(avg_power_taps);
 idx=find(avg_power_taps>strongest_tap-10); %keep taps with avg power > strongest tap -10dB

 ht=ht(idx,:);
 
 % fitting
 x=0:0.0001:2;
 y=zeros(size(x,2),size(idx,1));
 for i=1:size(idx,1)
    pd = fitdist(abs(fft(ht(i,:)))','Rician'); %en f ou t la rice ????
    y(:,i)=pdf(pd,x);
 end
 
 Legend=cell(1,size(idx,1));
 figure;
 for i=1:size(idx,1)
    plot(x,y(:,i));
    hold on;
    Legend{i}=strcat('tap', num2str(i));
 end
 grid on;
 ylabel('PDF');
 xlabel('|h|');
 legend(Legend);
 title('Rice fitting on channel impulse response taps LOS antenna 4');
 
  %%  Rice distrib fitting to channel impulse response taps NLOS
  
 ht_nlos = ifft(H_tensor_nlos(:,:,4));
 %ht_nlos=ht_nlos(1:256,:);
 
 htt_nlos=(abs(ht_nlos)).^2;
 avg_power_taps_nlos = 10*log10(mean(htt_nlos,2));
 strongest_tap_nlos = max(avg_power_taps_nlos);
 idx_nlos=find(avg_power_taps_nlos>strongest_tap_nlos-10);
 
 ht_nlos=ht_nlos(idx_nlos,:);
 
 y=zeros(size(x,2),size(idx_nlos,1));
 for i=1:size(idx_nlos,1)
    pd = fitdist(abs(fft(ht_nlos(i,:)))','Rician');
    y(:,i)=pdf(pd,x);
 end
 
 figure;
 for i=1:size(idx_nlos,1)
    plot(x,y(:,i));
    hold on;
    Legend{i}=strcat('tap', num2str(i));
 end
 grid on;
 ylabel('PDF');
 xlabel('|h|');
 legend(Legend);
 title('Rice fitting on channel impulse response taps NLOS antenna 4');
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average CIRs and PDP LOS

PDP = abs(autocorr(mean(ht,2))); %def of PDP

delta_tau=1/(2*params.B); %sampling time (tap duration)
tau=0:delta_tau:(length(PDP)-1)*delta_tau;
fun = mean((tau.^2).*PDP,'all');
taum =  mean((tau).*PDP,'all');
    
delay_spread = (sqrt(fun-taum^2)); %def of delay spread

Pmodel =PDP(1)*exp(-tau/(delay_spread)); %PDP modeled as exp
fun = mean((tau.^2).*Pmodel,'all');
taum =  mean((tau).*Pmodel,'all');
delay_spread_model = (sqrt(fun-taum^2));
coherence_BW = 1/(2*pi*delay_spread);

PDPverif=abs(ifft(autocorr(fft(mean(ht,2))))); %verif with autocorrelation in f domain
fun = mean((tau.^2).*PDPverif,'all');
taum =  mean((tau).*PDPverif,'all');
delay_spread_verif = (sqrt(fun-taum^2));
coherence_BW_verif = 1/(2*pi*delay_spread_verif);

figure;
plot(tau,Pmodel);
hold on;
plot(tau,PDP);
hold on;
grid on;
ylabel('PDP')
xlabel('\tau')
legend(strcat('PDP model: \sigma_T [\mus]= ',num2str(delay_spread_model*10e6)),strcat('PDP: \sigma_T [\mus]= ',num2str(delay_spread*10e6)))
title(strcat('Power Delay Profile LOS with  \sigma_{T_{verif}} [\mus]=', num2str(delay_spread_verif*10e6)))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average CIRs and PDP NLOS

PDPnlos = abs(autocorr(mean(ht_nlos,2)));
tau=0:delta_tau:(length(PDPnlos)-1)*delta_tau;

fun = mean((tau.^2).*PDPnlos,'all');
taum =  mean((tau).*PDPnlos,'all');
    
delay_spread_nlos = (sqrt(fun-taum^2));

Pmodel_nlos =PDPnlos(1)*exp(-tau/(delay_spread_nlos));
fun = mean((tau.^2).*Pmodel_nlos,'all');
taum =  mean((tau).*Pmodel_nlos,'all');
delay_spread_nlos_model = (sqrt(fun-taum^2));
coherence_BW_nlos = 1/(2*pi*delay_spread_nlos);

PDPverif_nlos=abs(ifft(autocorr(fft(mean(ht_nlos,2)))));
fun = mean((tau.^2).*PDPverif_nlos,'all');
taum =  mean((tau).*PDPverif_nlos,'all');
delay_spread_verif_nlos = (sqrt(fun-taum^2));
coherence_BW_verif_nlos = 1/(2*pi*delay_spread_verif_nlos);

figure;
plot(tau,Pmodel_nlos);
hold on;
plot(tau,PDPnlos);
hold on;
grid on;
ylabel('PDP')
xlabel('\tau')
legend(strcat('PDP model: \sigma_T [\mus]= ',num2str(delay_spread_nlos_model*10e6)),strcat('PDP: \sigma_T [\mus]= ',num2str(delay_spread_nlos*10e6)))
title(strcat('Power Delay Profile NLOS with \sigma_{T_{verif}} [\mus]=', num2str(delay_spread_verif_nlos*10e6)))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Influence of BW on PDP LOS

l=200; % length of PDP model
it=400; % Number of tests round(1638/it)
k=0;
figure;
for i=1:it:params.nActiveQ % Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor(1:i,:,4),params.nActiveQ);
    narrow_h=reshape(narrow_h,params.nActiveQ,[]);
    htt=(abs(narrow_h)).^2;
    avg_power_taps = 10*log10(sum(htt,2));
    strongest_tap = max(avg_power_taps);
    idx=find(avg_power_taps>strongest_tap-10);

    narrow_h=mean(narrow_h,2);
    narrow_h=narrow_h(idx,:);
    
    PDP_narrow = abs(autocorr(narrow_h));
    
    tau=0:delta_tau:(length(idx)-1)*delta_tau;
    fun = mean((tau.^2).*PDP_narrow,'all');
    taum =  mean((tau).*PDP_narrow,'all');
    
    delay_spread_narrow = (sqrt(abs(fun-taum^2))); %abs ? sinon 1 subcrr bug

    tau=0:delta_tau:(l-1)*delta_tau;
    Pmodel_narrow =PDP_narrow(1)*exp(-tau/(delay_spread_narrow));
 
    plot(tau,Pmodel_narrow);
    Legend{k}=strcat(num2str(i),' subcarriers');
    hold on;
end
grid on;
legend(Legend);
xlabel('\tau');
ylabel('PDP');
title('PDP for varying BW LOS')

k=0;
figure;
for i=1:it:params.nActiveQ % Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor(1:i,:,4),params.nActiveQ);
    narrow_h=reshape(narrow_h,params.nActiveQ,[]);
    %htt=(abs(narrow_h)).^2;
    %avg_power_taps = 10*log10(sum(htt,2));
    %strongest_tap = max(avg_power_taps);
    %idx=find(avg_power_taps>strongest_tap-10);

    %narrow_h=narrow_h(idx,:);
    
    tau=0:delta_tau:(params.nActiveQ-1)*delta_tau;
    plot(abs(mean(narrow_h,2)));
    hold on;
    Legend{k}=strcat(num2str(i),' subcarriers');
    title(sprintf('CIR'));
end
grid on;
legend(Legend);
xlabel('\tau');
ylabel('|h|');
title('CIR for varying BW LOS')

%% Influence of BW on PDP NLOS

k=0;
figure;
for i=1:it:params.nActiveQ % Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor_nlos(1:i,:,4),params.nActiveQ);
    narrow_h=reshape(narrow_h,params.nActiveQ,[]);
    htt=(abs(narrow_h)).^2;
    avg_power_taps = 10*log10(sum(htt,2));
    strongest_tap = max(avg_power_taps);
    idx=find(avg_power_taps>strongest_tap-10);

    narrow_h=mean(narrow_h,2);
    narrow_h=narrow_h(idx,:);
    
    PDP_narrow = abs(autocorr(narrow_h));
    
    tau=0:delta_tau:(length(idx)-1)*delta_tau;
    fun = mean((tau.^2).*PDP_narrow,'all');
    taum =  mean((tau).*PDP_narrow,'all');
    
    delay_spread_narrow = (sqrt(abs(fun-taum^2))); %abs ? sinon 1 subcrr bug

    tau=0:delta_tau:(l-1)*delta_tau;
    Pmodel_narrow =PDP_narrow(1)*exp(-tau/(delay_spread_narrow));
 
    plot(tau,Pmodel_narrow);
    Legend{k}=strcat(num2str(i),' subcarriers');
    hold on;
end
grid on;
legend(Legend);
xlabel('\tau');
ylabel('PDP');
title('PDP for varying BW NLOS')

k=0;
figure;
for i=1:it:params.nActiveQ % Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor_nlos(1:i,:,4),params.nActiveQ);
    narrow_h=reshape(narrow_h,params.nActiveQ,[]);
    %htt=(abs(narrow_h)).^2;
    %avg_power_taps = 10*log10(sum(htt,2));
    %strongest_tap = max(avg_power_taps);
    %idx=find(avg_power_taps>strongest_tap-10);

    %narrow_h=narrow_h(idx,:);
    
    tau=0:delta_tau:(params.nActiveQ-1)*delta_tau;
    plot(abs(mean(narrow_h,2)));
    hold on;
    Legend{k}=strcat(num2str(i),' subcarriers');
    title(sprintf('CIR'));
end
grid on;
legend(Legend);
xlabel('\tau');
ylabel('|h|');
title('CIR for varying BW NLOS')
