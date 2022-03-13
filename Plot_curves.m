clc
clear all
close all
% Declaration of SNR range, by rule of thumb, [-5dB,25dB] is enough.
snr=[-5:0.5:25];
x = zeros(length(snr),10);
y = size(x);
%% QAM Calculations
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',2);
    end
end
ber_qam_4 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4);
    end
end
ber_qam_16 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',6);
    end
end
ber_qam_64 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',8);
    end
end
ber_qam_256 = (mean(x,2)/480000).';

%%  PAM Calculations
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',1);
    end
end
ber_pam_2 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',2);
    end
end
ber_pam_4 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',3);
    end
end
ber_pam_8 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',4);
    end
end
ber_pam_16 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',5);
    end
end
ber_pam_32 = (mean(x,2)/480000).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',6);
    end
end
ber_pam_64 = (mean(x,2)/480000).';


for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',8);
    end
end
ber_pam_256 = (mean(x,2)/480000).';
%% Plotting PAMs
semilogy(snr,ber_qam_4,'-',snr,ber_qam_16,'-',snr,ber_qam_64,'-'); ylim([1e-5;1E0]);
figure;
semilogy(snr,ber_pam_2,'-');
hold on
semilogy(snr,ber_pam_4,'-');
semilogy(snr,ber_pam_8,'-');
semilogy(snr,ber_pam_16,'-');
semilogy(snr,ber_pam_32,'-');
semilogy(snr,ber_pam_64,'-');
semilogy(snr,ber_pam_256,'-');
ylim([1e-5 1]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('Bit Error Rate vs E_B/N_0 (dB) for Different Modulation Types')
grid on;
legend('2-PAM','4-PAM','8-PAM','16-PAM','32-PAM','64-PAM','256-PAM')
% 
% 
%% Plotting QAMs
%semilogy(snr,ber_qam_4,'-',snr,ber_qam_16,'-',snr,ber_qam_64,'-'); ylim([1e-5;1E0]);
figure;
semilogy(snr,ber_qam_4,'-');
hold on
semilogy(snr,ber_qam_16,'-');
semilogy(snr,ber_qam_64,'-');
semilogy(snr,ber_qam_256,'-');
ylim([1e-5 1]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('Bit Error Rate vs E_B/N_0 (dB) for Different Modulation Types')
grid on;
legend('4-QAM','16-QAM','64-QAM','256-QAM')
%% Illustration of Nyquist Filter No ISI
filter_taps_t = conv(filter_taps_t,filter_taps_t);
figure;
yaxis =  (-(length(filter_taps_t)-1)/2)*1/8e6 : 1/8e6 : ((length(filter_taps_t)-1)/2)*1/8e6;
p = plot(yaxis,filter_taps_t,'g--o');
hold on
plot(yaxis,filter_taps_t);
title("Nyquist Filter (RRC) Impulse Response")
p.MarkerSize = 8;
p.MarkerFaceColor = "black";
p.Color = "green";
legend("Impulse Response Sampled at Symbol Durations","Full Impulse Response")
p.MarkerIndices = 1:4:length(filter_taps_t);
ylabel("Normalized Amplitude (such that 1/2 RRC Filter has unit energy)")
xlabel("Time (seconds)")
grid on
%% Illustration of Nyquist Filter is Bandlimited

% 1. method: Using only FFT method
% figure;
% F_samp = 8e6;
% plot(linspace(-F_samp/2,F_samp/2,length(filter_taps_f)),filter_taps_f)
% title("Fourier Transform of Nyquist Filter")
% xlabel("Frquency (Hz)")

% 2. method: Using FFT and making PSD
% Fs = F_samp;
% t = 0:1/Fs:1-1/Fs;
%x = filter_taps_t;
%y = upsampled_symbols;
% N = length(x);
% xdft = fft(x);
% xdft = xdft(1:N/2+1);
% psdx = (1/(Fs*N)) * abs(xdft).^2;
% psdx(2:end-1) = 2*psdx(2:end-1);
% freq = 0:Fs/length(x):Fs/2;
% 
% plot(freq,10*log10(psdx))
% grid on
% title('Periodogram Using FFT')
% xlabel('Frequency (Hz)')
% ylabel('Power/Frequency (dB/Hz)')

% 3. Method: Using Periodogram function
figure;
Fs = 8e6;
x = filter_taps_t;
periodogram(x,rectwin(length(x)),length(x),Fs,"centered")
title("PSD Estimation of filter")

figure;
y = upsampled_symbols;
periodogram(y,rectwin(length(y)),length(y),Fs,"centered")
title("PSD Estimation of Non-Filtered Upsampled Symbols")

figure;
z = filtered_symbols;
periodogram(z,rectwin(length(z)),length(z),Fs,"centered")
title("PSD Estimation of Filtered Symbols")
