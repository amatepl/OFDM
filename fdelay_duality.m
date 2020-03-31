clc; close all; clear;
addpath("MA2_libs")
params = load('TestParam.mat').TestParam;

HLOS = load("H_LOS.mat").H;
HNLOS = load("H_NLOS.mat").H;

H_tensor=zeros(params.nActiveQ,99,4);
H_tensor_nlos=zeros(params.nActiveQ,99,4);

for i = 1:4
    H_LOS = fftshift(HLOS(i,:,:));
    H_NLOS = fftshift(HNLOS(i,:,:));
    for j = 1:99
        H_tensor(:,j,i) = H_LOS(1,params.ActiveQIndex,j);
        H_tensor_nlos(:,j,i) = H_NLOS(1,params.ActiveQIndex,j);
    end
end

 
 z(:,:)=H_tensor(:,:,2);
 w(:,:)=H_tensor_nlos(:,:,2);
 figure;
 subplot(1,2,1)
 plot(abs(z))
 grid on;
 ylabel('|H|')
 xlabel('subcarriers')
 title('LOS')
 subplot(1,2,2)
 plot(abs(w))
 grid on;
 ylabel('|H|')
 xlabel('subcarriers')
 title('NLOS')
 sgtitle('Channel frequency response estimates')
%  
%  v(:,:)=H_tensor(:,10,:);
%  u(:,:)=H_tensor_nlos(:,10,:);
%  figure;
%  subplot(1,2,1)
%  plot(abs(v))
%  grid on;
%  ylabel('|h|')
%  xlabel('Active subcarriers')
%  title('LOS')
%  subplot(1,2,2)
%  plot(abs(u))
%  grid on;
%  ylabel('|h|')
%  xlabel('Active subcarriers')
%  title('NLOS')
%  sgtitle('Channel frequency response estimates')
%  
%  x_val=0:0.001:0.5;
%  a(:,1)=sum(H_tensor(:,10,:),1);
%  b(:,1)=sum(H_tensor_nlos(:,10,:),1);
%  p=fitdist(abs(a),'Rician');
%  q=fitdist(abs(b),'Rician');
%  y=pdf(p,x_val);
%  w=pdf(q,x_val);
%  K=10*log10(p.s^2/2*p.sigma^2); %rice factor
%  tot= p.s^2 + 2*p.sigma^2;
%  Kq=10*log10(q.s^2/2*q.sigma^2); %rice factor
%  figure;
%  plot(x_val,y)
%  figure;
%  plot(x_val,w)
 
 %% fit rice distrib to channel frequency response for chosen antenna and chosen frame
 x_val=0:0.01:5; %values of |h| where the pdf is evaluated
 ylos=zeros(size(x_val,2),4); ynlos=zeros(size(x_val,2),4);
 Klos=zeros(4); Knlos=zeros(4);
 for i=1:4
    a(:,:)=H_tensor(:,10,i);
    b(:,:)=H_tensor_nlos(:,10,i);
    p1 = fitdist(abs(a),'Rician');
    p2 = fitdist(abs(b),'Rician');
    ylos(:,i)=pdf(p1,x_val);
    ynlos(:,i)=pdf(p2,x_val);
    Klos(i)=10*log10(p1.s^2/2*p1.sigma^2); %rice factor
    Knlos(i)=10*log10(p2.s^2/2*p1.sigma^2);
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
 %%  Rice distrib fitting to channel impulse response taps NLOS
 % taps distributions of CIR's
 
 ht2 = ifft((H_tensor_nlos(:,:,2)));
 ht2 = ifftshift(ht2);
 ht2=ht2(820:end,:);

 %alignement
  [~,I2] = max(abs(ht2).^2);
  httn2=zeros(819,length(I2));
  for i=1:length(I2)
      httn2(1:end-I2(i)+1,i)=ht2(I2(i):end,i);
  end
  
 figure;
 plot(abs(httn2));
 ylabel('|h|')
 title('Channel realizations CIRs NLOS');
 
 
  PDP2=zeros(819,99);
  for i=1:99
     PDP2(:,i) = abs(autocorr(httn2(:,i)));
  end
 PDP2=mean(PDP2,2);
%PDP=mean(abs(httn).^2,2);
delta_tau=1/(2*params.B); 
tau=0:delta_tau:delta_tau*(length(PDP2)-1);

figure;
plot(tau,(PDP2));
grid on;
xlabel('\tau')
ylabel('PDP [dB]')
title('PDP over all CIR taps NLOS')

avg_power_taps2=10*log10(PDP2);
strongest_tap2=10*log10(PDP2(1));
idx2=find(avg_power_taps2>strongest_tap2-10);
 
PDP2=PDP2(idx2);
htt2=abs(httn2(idx2,:));
 
 % fitting
 x=0:0.0001:3;
 y2=zeros(size(x,2),length(idx2));
 K_taps2=zeros(length(idx2));
 for i=1:length(idx2)
    pd2 = fitdist((htt2(i,:)'),'Rician');
    K_taps2(i)=10*log10(pd2.s^2/2*pd2.sigma^2);
    y2(:,i)=pdf(pd2,x);
 end
 
 Legend=cell(1,length(idx2));
 figure;
 for i=1:length(idx2)
    plot(x,y2(:,i));
    hold on;
    Legend{i}=strcat('tap', num2str(i), ' KdB= ', num2str(K_taps2(i)));
 end
 grid on;
 ylabel('PDF');
 xlabel('|h|');
 legend(Legend);
 title('Rice fitting on channel impulse response taps NLOS antenna 2');
 
 %% PDP - delay spread
 tau=0:delta_tau:delta_tau*(length(PDP2)-1);
% 
P2=sum(PDP2);
tum=sum(tau'.*PDP2)/P2;
fum=sum((tau.^2)'.*PDP2)/P2;

delay_spread2=sqrt(fum-tum^2);
coherence_BW2 = 1/(2*pi*delay_spread2);

Pmodel2(:,1) =PDP2(1)*exp(-tau/(delay_spread2)); %PDP modeled as exp
P2=sum(Pmodel2);
tum=sum(tau'.*Pmodel2)/P2;
fum=sum((tau.^2)'.*Pmodel2)/P2;
delay_spread_model2=sqrt(fum-tum^2);
coherence_BW_model2 = 1/(2*pi*delay_spread_model2);

R2=zeros(1638,99);
for i=1:99
    R2(:,i)=autocorr((H_tensor_nlos(:,i,2)));
end
R2=mean(R2,2);

figure;
plot(abs(R2));
title('Channel Frequency correlation NLOS')
ylabel('|R(\Delta f)|')
xlabel('\Delta f')

PDPverif2=(ifft(R2)); %verif with autocorrelation in f domain

PDPverif2=ifftshift(PDPverif2);
PDPverif2=abs(PDPverif2(820:end));
[~,id2]=max(PDPverif2);
PDPverif2=PDPverif2(id2:end);
PDPverif2=PDPverif2(idx2);

tau=0:delta_tau:delta_tau*(length(PDPverif2)-1);

P2=sum(PDPverif2);
tum=sum(tau'.*PDPverif2)/P2;
fum=sum((tau.^2)'.*PDPverif2)/P2;

delay_spread_verif2=sqrt(fum-tum^2);
coherence_BW_verif2 = 1/(2*pi*delay_spread_verif2);

tau=0:delta_tau:delta_tau*(length(PDP2)-1);
figure;
plot(tau,Pmodel2);
hold on;
plot(tau,PDP2);
hold on;
grid on;
ylabel('PDP')
xlabel('\tau')
legend(strcat('PDP model:  \sigma_T [\mus]= ',num2str(delay_spread_model2*1e6),' and \Delta f_c [MHz]= ',num2str(coherence_BW_model2*1e-6)),strcat('PDP: \sigma_T [\mus]= ',num2str(delay_spread2*1e6),'and \Delta f_c [MHz]= ',num2str(coherence_BW2*1e-6)))
title(strcat('Power Delay Profile NLOS with  \sigma_{T_{verif}} [\mus]=', num2str(delay_spread_verif2*1e6)))

%%  Rice distrib fitting to channel impulse response taps NLOS
 % taps distributions of CIR's
 
%  figure;
%  plot(abs(H_tensor_nlos(:,1,2)));
 
 ht = ifft((H_tensor(:,:,2)));
 ht = ifftshift(ht);
 ht=ht(820:end,:);


 %alignement
  [~,I] = max(abs(ht).^2);
  httn=zeros(819,length(I));
  for i=1:length(I)
      httn(1:end-I(i)+1,i)=ht(I(i):end,i);
  end
  
 figure;
 plot(abs(httn));
 ylabel('|h|')
 title('Channel realizations CIRs LOS');
 
 
 PDP=zeros(819,99);
 for i=1:99
    PDP(:,i) = abs(autocorr(httn(:,i)));
 end
PDP=mean(PDP,2);
%PDP=mean(abs(httn).^2,2);

delta_tau=1/(params.B); 
tau=0:delta_tau:delta_tau*(length(PDP)-1);

figure;
plot(tau,(PDP));
grid on;
xlabel('\tau')
ylabel('PDP [dB]')
title('PDP over all CIR taps LOS')

avg_power_taps=10*log10(PDP);
strongest_tap=10*log10(max(PDP));
idx=find(avg_power_taps>strongest_tap-10);
 
PDP=PDP(idx);
htt=abs(httn(idx,:));
 
 % fitting
 x=0:0.0001:2;
 y=zeros(size(x,2),length(idx));
 K_taps=zeros(length(idx));
 for i=1:length(idx)
    pd = fitdist((htt(i,:)'),'Rician');
    K_taps(i)=10*log10(pd.s^2/2*pd.sigma^2);
    y(:,i)=pdf(pd,x);
 end
 
 Legend=cell(1,length(idx));
 figure;
 for i=1:length(idx)
    plot(x,y(:,i));
    hold on;
    Legend{i}=strcat('tap', num2str(i), ' KdB= ', num2str(K_taps(i)));
 end
 grid on;
 ylabel('PDF');
 xlabel('|h|');
 legend(Legend);
 title('Rice fitting on channel impulse response taps LOS antenna 2');
 
 %% PDP - delay spread NLOS
tau=0:delta_tau:delta_tau*(length(PDP)-1);

P=sum(PDP);
tum=sum(tau'.*PDP)/P;
fum=sum((tau.^2)'.*PDP)/P;

delay_spread=sqrt(fum-tum^2);
coherence_BW = 1/(2*pi*delay_spread);

Pmodel(:,1) =PDP(1)*exp(-tau/(delay_spread)); %PDP modeled as exp
P=sum(Pmodel);
tum=sum(tau'.*Pmodel)/P;
fum=sum((tau.^2)'.*Pmodel)/P;
delay_spread_model=sqrt(fum-tum^2);
coherence_BW_model = 1/(2*pi*delay_spread_model);

R=zeros(1638,99);
for i=1:99
    R(:,i)=autocorr((H_tensor(:,i,2)));
end
R=mean(R,2);

figure;
plot(abs(R));
title('Channel Frequency correlation LOS')
ylabel('|R(\Delta f)|')
xlabel('\Delta f')

PDPverif=(ifft(R)); %verif with autocorrelation in f domain

PDPverif=ifftshift(PDPverif);
PDPverif=abs(PDPverif(820:end));
%[~,id]=max(PDPverif);
%PDPverif=PDPverif(id:end);
PDPverif=PDPverif(idx);

tau=0:delta_tau:delta_tau*(length(PDPverif)-1);

P=sum(PDPverif);
tum=sum(tau'.*PDPverif)/P;
fum=sum((tau.^2)'.*PDPverif)/P;

delay_spread_verif=sqrt(fum-tum^2);
coherence_BW_verif = 1/(2*pi*delay_spread_verif);

tau=0:delta_tau:delta_tau*(length(PDP)-1);
figure;
plot(tau,Pmodel);
hold on;
plot(tau,PDP);
hold on;
grid on;
ylabel('PDP')
xlabel('\tau')
legend(strcat('PDP model:  \sigma_T [\mus]= ',num2str(delay_spread_model*1e6),' and \Delta f_c [MHz]= ',num2str(coherence_BW_model*1e-6)),strcat('PDP: \sigma_T [\mus]= ',num2str(delay_spread*1e6),' and \Delta f_c [MHz]= ',num2str(coherence_BW*1e-6)))
title(strcat('Power Delay Profile LOS with  \sigma_{T_{verif}} [\mus]=', num2str(delay_spread_verif*1e6)))
 
%% Influence of BW on PDP NLOS

k=0;
it=[1,40,400,900,1638];
l=100;
figure;
for i=it% Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor_nlos(1:i,:,2),[],1); 
    narrow_h = ifftshift(narrow_h);
    narrow_h=narrow_h(floor(size(narrow_h,1)/2)+1:end,:);
    %narrow_h=[narrow_h;zeros(819-i,99)];

    
    %alignement
    [~,id] = max(abs(narrow_h).^2,[],1);
    httn=zeros(size(narrow_h,1),length(id));
    for j=1:length(id)
        httn(1:end-id(j)+1,j)=narrow_h(id(j):end,j);
    end
    
    %figure;
    %stem(mean(abs(httn),2))

     PDP_narrow=zeros(size(narrow_h,1),99);
    for j=1:99
     PDP_narrow(:,j) = abs(autocorr(httn(:,j)));
    end
    PDP_narrow=mean(PDP_narrow,2);
    
    avg_power_taps = 10*log10(PDP_narrow);
    strongest_tap = 10*log10(max(PDP_narrow));
    idx=find(avg_power_taps>strongest_tap-10);

    PDP_narrow=PDP_narrow(idx);
    
    %delta_tau=1/((params.B/1638)*i);
    tau=0:delta_tau:delta_tau*(length(PDP_narrow)-1);

    P_narrow=sum(PDP_narrow*delta_tau);
    tum=sum(tau'.*PDP_narrow*delta_tau)/P_narrow;
    fum=sum((tau.^2)'.*PDP_narrow*delta_tau)/P_narrow;
    
    delay_spread_narrow=sqrt(fum-tum^2);
    coherence_BW_narrow = 1/(2*pi*delay_spread_narrow);
    
    Pmodel_narrow =PDP_narrow(1)*exp(-tau/(delay_spread_narrow));
 
    
    plot(tau,Pmodel_narrow);
    hold on;
    Legend{k}=strcat(num2str(i),' subcarriers with \sigma_T [\mus]',num2str(delay_spread_narrow*1e6));
end
grid on;
legend(Legend);
xlabel('\tau');
ylabel('PDP');
title('PDP for varying BW NLOS')

k=0;
it=[1,10,800,1638];
figure;
for i=it% Number of subcarrier
    k=k+1;
    narrow_h = ifft(H_tensor_nlos(1:i,:,2),[],1); 
    narrow_h = ifftshift(narrow_h);
    narrow_h=narrow_h(floor(size(narrow_h,1)/2)+1:end,:);
    narrow_h=[narrow_h;zeros(819-i,99)];
    
    %alignement
    [~,id] = max(abs(narrow_h).^2,[],1);
    httn=zeros(size(narrow_h,1),length(id));
    for j=1:length(id)
        httn(1:end-id(j)+1,j)=narrow_h(id(j):end,j);
    end

    h_narrow=mean(abs(httn),2);
    tau=0:delta_tau:delta_tau*(length(h_narrow)-1);
 
    subplot(2,2,k);
    stem(tau,h_narrow);
    hold on;
    legend(strcat(num2str(i),' subcarriers'));
    grid on;
    xlabel('\tau');
    ylabel('PDP');
end
sgtitle('CIRs for varying BW NLOS')