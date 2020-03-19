function y = testFunc(signal_tx,H)

H = fftshift(H(1,:)); % Take only one impulse respnse

h = ifft((H),[],2); % Impulse response in time domain

% h(1:512) = h(1:512) + flip(h(end - 512 +1:end));

%% WORKING IN TIME DOMAIN

% Convolution of the signal with the impulse response
% signal_conv = conv(signal_tx(1:2048+512),h(1:100),'full');
N = 2048;
o = zeros(1,2048);
o(1:512) = h(1:512);
signal_conv = cconv(h,signal_tx(1,513:2048+512),N);
signal_conv = conv(h(1:512),signal_tx(1,1:2048+512));

%% WORKING IN FREQUENCY DOMAIN

% CP removal
signal = reshape(signal_tx,[],32);
signal = signal(513:end,:);

% Frequency domain
Signal = fft(signal,[],1);
H = ifft(H);
H = fft(H(1:512),2048);
Signal = Signal.*(H.');

% Back to time domain
signal = ifft(Signal,[],1);

%% PLOTS

% Plot h
figure, hold on;
plot(abs(h));
xlabel('Time');
title('h');

% Plot of the signal_tx
figure, hold on;
plot(abs(signal_tx(513:2048+512)));
xlabel('Time');
title('Signal tx');

% Plot of the signal convolved with the channel in time domain
m = 0;

figure, hold on;
subplot(2,1,1)
plot(abs(signal_conv(513+m:(2048+512)*1+m)));
% plot(abs(signal_conv(:)));
xlabel('Time');
title('Signal convolved in time domain');

subplot(2,1,2)
plot(abs(fft(signal_conv(513+m:(2048+512)*1+m))));
% plot(abs(fft(signal_conv,[],2)));
xlabel('Frequecies');
title('Signal convolved in time domain');

% Plot of the signal multiplied with the channel in frequency domain
figure, hold on;
subplot(2,1,1)
plot(abs(signal(:,1)));
xlabel('Time');
title('Signal multiplied in frequency domain');

subplot(2,1,2)
plot(abs(Signal(:,1)));
xlabel('Frequecies');
title('Signal multiplied in frequency domain');


keyboard