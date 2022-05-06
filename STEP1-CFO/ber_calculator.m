function [error_symbols] = ber_calculator(snr,modulation_type,Nbps,ppm)

number_of_bits = 1e6;
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);
index =2000;
desired_snr = 15;
desired_ppm = 10;
fc = 2e+9; % carrier frquency

cfo = ppm*1e-6*fc;

f_cut = 1e6; % cutoff frquency
F_sampling = 8e6; %sampling frequency, it has been increase significantly to 
% be able to simulate time shift
T_sampling = 1/F_sampling;
F_symbol =2*f_cut; %symbol frequency
T_symbol = 1/F_symbol;
upsample_rate = F_sampling/F_symbol;
rand_phase = 1/8; %Will be multiplied with 2*pi afterwards
upsampled_symbols = upsample(transmitted_symbols,upsample_rate);
if snr == desired_snr && ppm == desired_ppm
figure;
a = sgtitle('Scatterplots of Different Samples Under SNR=15 dB, CFO = 10 ppm, Phase Shift=\pi/4 rads (Samples are downsampled with factor 500)','Color','black');
a.FontSize = 16;
a.FontWeight = 'bold';

subplot(2,2,1)
scatter(real(upsampled_symbols),imag(upsampled_symbols),5,'filled','b')
title("QAM Symbols")
end
beta=0.8;
span = 2;
taps = upsample_rate*span+1;

% Applying first filter
[filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);


E_symbol_av = sum(abs(transmitted_symbols.*transmitted_symbols))/length(transmitted_symbols);

%Energy Per Bit, To find it, I divide symbol energy into bits per symbol
Eb = E_symbol_av/Nbps;

% Adding noise to symbols
noisy_symbols = add_noise(snr,modulation_type,Eb,filtered_symbols);
if snr == desired_snr && ppm == desired_ppm
subplot(2,2,2)
scatter(real(noisy_symbols(1:index:end)),imag(noisy_symbols(1:index:end)),5,'filled','b')
title("Noise Added QAM Symbols")
end
%Adding CFO
t1 = ((0:length(noisy_symbols)-1))*(1/(F_sampling));

for i = 1:length(t1)
    noisy_symbols(i) = noisy_symbols(i)*exp(1j*(2*pi* (cfo*t1(i)+rand_phase) ));
end
if snr == desired_snr && ppm == desired_ppm
subplot(2,2,3)
scatter(real(noisy_symbols(1:index:end)),imag(noisy_symbols(1:index:end)),5,'filled','b')
title("CFO Added QAM Symbols")
end
%  Applying second filter
filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

% Eliminating filtered data
filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);



%Compensate CFO
t2 = ((0:length(filtered_noisy_symbols)-1))*(1/(F_sampling));
 
for i = 1:length(t2)
    filtered_noisy_symbols(i) = filtered_noisy_symbols(i)*exp(-1j*(2*pi* (cfo*t1(i)+rand_phase) ));
end
if snr == desired_snr && ppm == desired_ppm
subplot(2,2,4)
scatter(real(filtered_noisy_symbols(1:index:end)),imag(filtered_noisy_symbols(1:index:end)),5,'filled','b')
title("CFO Compensated QAM Symbols")
end

% Downsampling
filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

% Demapping
demodulated_bits = demapping(filtered_noisy_symbols,Nbps,modulation_type);


error_symbols = sum(bits ~= demodulated_bits);


filter_tapss = filter_taps_in_time;
filter_tapss_freq=filter_taps_in_freq;
end


