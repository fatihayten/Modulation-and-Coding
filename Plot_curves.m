clc
clear all
close all

% Declaration of SNR range, by rule of thumb, [-5dB,25dB] is enough.
snr=(-5:0.5:25);

% For each SNR value, I will simulate 10 times, then I will calculate the
% average BER, because only 1 simulation may result in failure. 

% If you dont trust on your RAM/Processor, you can decrease it to 3 or 4.

% I will store the BER values in an empty matrix.
x = zeros(length(snr),10); 

% I am calculating the size of x, to do for loop operation.
y = size(x); 
%% QAM Calculations

% In each for loop, I am calculating BER for different modulation schemes. 

% 4-QAM
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',2);
    end
end
% To find mean of BER for different simulations, I am calculating the mean.
ber_qam_4 = (mean(x,2)/480000).';

% 16-QAM
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4);
    end
end
ber_qam_16 = (mean(x,2)/480000).';

% 64-QAM
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',6);
    end
end
ber_qam_64 = (mean(x,2)/480000).';

% 256-QAM
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',8);
    end
end
ber_qam_256 = (mean(x,2)/480000).';
%% Theoretical BER Curves of QAM Using Embedded MATLAB Functions

 
% In the project, we supposed to compare the experimental results with
% theoretical ones. To do so, I found some functions on MATLAB website.
% Reference: https://nl.mathworks.com/help/comm/ref/biterr.html


% Defining the Eb/No range, should be same as the vector in Line 6.
EbNoVec = (-5:0.5:25)'; 
ber_qam_4_theoretical = berawgn(EbNoVec,'qam',4); % 4-QAM BER Curve
ber_qam_16_theoretical = berawgn(EbNoVec,'qam',16); % 4-QAM BER Curve
ber_qam_64_theoretical = berawgn(EbNoVec,'qam',64); % 4-QAM BER Curve
ber_qam_256_theoretical = berawgn(EbNoVec,'qam',256); % 4-QAM BER Curve

%%  PAM Calculations
% IN PROJECT, WE ONLY CONSIDER QAM MODULATION


% % I am doing the same steps for PAM modulation.
% 
% % 2-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',1);
%     end
% end
% ber_pam_2 = (mean(x,2)/480000).';
% 
% % 4-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',2);
%     end
% end
% ber_pam_4 = (mean(x,2)/480000).';
% 
% % 8-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',3);
%     end
% end
% ber_pam_8 = (mean(x,2)/480000).';
% 
% % 16-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',4);
%     end
% end
% ber_pam_16 = (mean(x,2)/480000).';
% 
% % 32-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',5);
%     end
% end
% ber_pam_32 = (mean(x,2)/480000).';
% 
% % 64-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',6);
%     end
% end
% ber_pam_64 = (mean(x,2)/480000).';
% 
% % 256-PAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'pam',8);
%     end
% end
% ber_pam_256 = (mean(x,2)/480000).';

%% Theoretical BER Curves of PAM Using Embedded MATLAB Functions
% IN PROJECT, WE ONLY CONSIDER QAM MODULATION


% EbNoVec = (-5:0.5:25)'; 
% ber_pam_2_theoretical = berawgn(EbNoVec,'pam',2); % 4-QAM BER Curve
% ber_pam_4_theoretical = berawgn(EbNoVec,'pam',4); % 4-QAM BER Curve
% ber_pam_8_theoretical = berawgn(EbNoVec,'pam',8); % 4-QAM BER Curve
% ber_pam_16_theoretical = berawgn(EbNoVec,'pam',16); % 4-QAM BER Curve
% ber_pam_32_theoretical = berawgn(EbNoVec,'pam',32); % 4-QAM BER Curve
% ber_pam_64_theoretical = berawgn(EbNoVec,'pam',64); % 4-QAM BER Curve
% ber_pam_256_theoretical = berawgn(EbNoVec,'pam',256); % 4-QAM BER Curve

%% Plotting PAM BER Curves
% IN PROJECT, WE ONLY CONSIDER QAM MODULATION


% % In plotting parts of PAM and QAM, by using "hold on" command, I am
% % plotting different modulation types and theoretical curves one on
% % the top of the other
% 
% figure;
% % Experimental results
% semilogy(snr,ber_pam_2,'-');
% hold on
% semilogy(snr,ber_pam_4,'-');
% semilogy(snr,ber_pam_8,'-');
% semilogy(snr,ber_pam_16,'-');
% semilogy(snr,ber_pam_32,'-');
% semilogy(snr,ber_pam_64,'-');
% semilogy(snr,ber_pam_256,'-');
% 
% % Theoretical results
% semilogy(EbNoVec,ber_pam_2_theoretical,'.')
% semilogy(EbNoVec,ber_pam_4_theoretical,'.')
% semilogy(EbNoVec,ber_pam_8_theoretical,'.')
% semilogy(EbNoVec,ber_pam_16_theoretical,'.')
% semilogy(EbNoVec,ber_pam_32_theoretical,'.')
% semilogy(EbNoVec,ber_pam_64_theoretical,'.')
% semilogy(EbNoVec,ber_pam_256_theoretical,'.')
% 
% % Indicating the BER range
% ylim([1e-5 1])
% xlabel('E_B/N_0 [dB]');
% ylabel('BER');
% title(['Comparison of Experimental and Theoretical Results of Bit Error...' ...
%     ' Rate vs E_B/N_0 (dB) for Different Modulation Types'])
% grid on;
% legend('2-PAM Experimental','4-PAM Experimental','8-PAM Experimental','16-PAM Experimental','32-PAM Experimental','64-PAM Experimental','256-PAM Experimental','2-PAM Theoretical','4-PAM Theoretical','8-PAM Theoretical','16-PAM Theoretical','32-PAM Theoretical','64-PAM Theoretical','256-PAM Theoretical')

%% Plotting QAM BER Curves
%semilogy(snr,ber_qam_4,'-',snr,ber_qam_16,'-',snr,ber_qam_64,'-'); ylim([1e-5;1E0]);

figure;
% Experimental results
semilogy(snr,ber_qam_4,'-');
hold on
semilogy(snr,ber_qam_16,'-');
semilogy(snr,ber_qam_64,'-');
semilogy(snr,ber_qam_256,'-');

% Theoretical results
semilogy(EbNoVec,ber_qam_4_theoretical,'.')
semilogy(EbNoVec,ber_qam_16_theoretical,'.')
semilogy(EbNoVec,ber_qam_64_theoretical,'.')
semilogy(EbNoVec,ber_qam_256_theoretical,'.')

% Indicating the BER range
ylim([1e-5 1]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title(['Comparison of Experimental and Theoretical Results of Bit Error',...
    ' Rate vs E_B/N_0 (dB) for Different Modulation Types']);
grid on;
legend('4-QAM Experimental','16-QAM Experimentaml','64-QAM Experimental',...
    '256-QAM Experimental','4-QAM Theoretical','16-QAM Theoretical',...
    '64-QAM Theoretical','256-QAM Theoretical')
%% Illustration of Nyquist Filter Has No ISI

% Convolving Half RRC with itself, so find the full RRC filter
filter_taps_t = conv(filter_taps_t,filter_taps_t);

% Now plot the result, since sps is 4, starting from the center point,
% at each 4 samples for right and left, the sampled will zero (it is 
% a theoretical assumption, in MATLAB it is about e-5 or e-6)

% In the plot, the samples that are indicated with dots corresponds to 
% time axis at nT (T=symbol duration=T_sampling*upsample_rate).
figure;
% Arranging the y axis properly to illustrate No ISI
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

% In this step, I found 3 different method, but the coolest one is third
% one, because it uses periodogram function and this function is super
% easy to use.

% You will see that, cut-off frequency of filter is almost 1 MHz.
% Also; non filtered symbols covers the whole frequency band which
% is too bad for real life implementation.
% The filtered symbols covers only 2 MHz bandwidth, just like the filter.
% So, symbols are filtered successfully.

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
Fs = 8e6; % Sampling Frequency
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


