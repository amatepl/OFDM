 clc; close all; clear;
 
 params = load('TestParam.mat').TestParam;

 H_tensor=load('H_tensor.mat').H_tensor;
 H_tensor_nlos=load('H_tensor_nlos.mat').H_tensor_nlos;

 x_val=0:0.01:5;
 ylos=zeros(size(x_val,2),4); ynlos=zeros(size(x_val,2),4);
 Klos=zeros(4); Knlos=zeros(4);
 for i=1:4
    p1 = fitdist(abs(H_tensor(:,10,i)),'Rician');
    p2 = fitdist(abs(H_tensor_nlos(:,10,i)),'Rician');
    ylos(:,i)=pdf(p1,x_val);
    ynlos(:,i)=pdf(p2,x_val);
    Klos(i)=10*log10(p1.s^2/2*p1.sigma);
    Knlos(i)=10*log10(p2.s^2/2*p1.sigma);
 end
 
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
 
 ht = ifft(H_tensor(:,:,4));
 %ht=ht(1:512,:);
 
 htt=(abs(ht)).^2;
 avg_power_taps = 10*log10(sum(htt,2));
 strongest_tap = max(avg_power_taps);
 idx=find(avg_power_taps>strongest_tap-10);

 ht=ht(idx,:);
 
 
 x=0:0.0001:2;
 y=zeros(size(x,2),size(idx,1));
 for i=1:size(idx,1)
    pd = fitdist(abs(fft(ht(i,:)))','Rician');
    y(:,i)=pdf(pd,x);
 end
 
 
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
 
 
 ht_nlos = ifft(H_tensor_nlos(:,:,4));
 %ht_nlos=ht_nlos(1:512,:);
 
 htt_nlos=(abs(ht_nlos)).^2;
 avg_power_taps_nlos = 10*log10(sum(htt_nlos,2));
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

ht = mean(ht,2);
PDP = abs(autocorr(ht));

delta_tau=1/(2*params.B);
tau=0:delta_tau:(length(PDP)-1)*delta_tau;
P= sum(PDP,'all');
fun = mean((tau.^2).*PDP,'all');
taum =  mean((tau).*PDP,'all');
    
delay_spread = (sqrt(fun-taum^2));

Pmodel =PDP(1)*exp(-tau/(delay_spread));
fun = mean((tau.^2).*Pmodel,'all');
taum =  mean((tau).*Pmodel,'all');
delay_spread_model = (sqrt(fun-taum^2));
coherence_BW = 1/(2*pi*delay_spread_model);

figure;
plot(Pmodel);
hold on;
plot(PDP);
hold on;
grid on;
ylabel('PDP')
xlabel('tau')
legend('PDP model','PDP')
title('Power Delay Profile')

R=autocorr(fft(ht));
PDPverif=abs(ifft(R));
fun = mean((tau.^2).*PDPverif,'all');
taum =  mean((tau).*PDPverif,'all');
delay_spread_verif = (sqrt(fun-taum^2));
coherence_BW_verif = 1/(2*pi*delay_spread_verif);
 