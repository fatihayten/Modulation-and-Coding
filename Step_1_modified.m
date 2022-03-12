clc
close all
number_of_bits = 6000;
bits = random_bit_generator(number_of_bits);
Nbps = 4;
modulation_type = "qam";
transmitted_symbols = mapping(bits,Nbps,modulation_type);
%scatterplot(transmitted_symbols,10)
%plot(transmitted_symbols,'g*')
F_sampling = 8e6; %sampling frequency
F_symbol = 2e6; %symbol frequency
upsample_rate = F_sampling/F_symbol; %Should be integer, because L argument of...
% upsample function should be integer-valued.
% https://nl.mathworks.com/help/control/ref/lti.upsample.html?s_tid=doc_ta

upsampled_symbols = upsample(transmitted_symbols,upsample_rate);

beta=0.3;
span = 10;
taps = upsample_rate*span+1;
a = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
figure;
plot(a)
hold on

%checking the difference between self-written filter and rcosdesign filter
b = rcosdesign(beta,span,upsample_rate,"normal");
plot(b)
legend("Self-written filter","rcosdesing filter")
difference_of_2_filters = sum(abs(a-b))
title("Nyquist Filter in Time Domain")

%illustration of cancellation of the inter-symbol interference
% figure;
% hold on
% plot(conv(a,a))
% plot(circshift(conv(a,a),upsample_rate))
% plot(circshift(conv(a,a),2*upsample_rate))

% Passing the symbols through the filter1
filtered_symbols = conv(upsampled_symbols,a);


% Average Symbol Energy
E_symbol_av = sum(abs(transmitted_symbols))/length(transmitted_symbols);

%Energy Per Bit
Eb = E_symbol_av/Nbps;

% Declaration of SNR range, by rule of thumb, [-5dB,25dB] is enough.
snr_range = -5:5;

% Adding noise to symbols
noisy_symbols = add_noise(snr_range,modulation_type,Eb,filtered_symbols,F_sampling);

% Passing the symbols through the filter1
filtered_noisy_symbols = conv(noisy_symbols,a);

% Eliminating filtered data
filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);

% Downsampling
filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

% Demapping
demodulated_bits = (demapping(filtered_noisy_symbols,Nbps,modulation));

% Finding The Probability of Error in Bits
number_of_different_bits = sum(bits ~= demodulated_bits)


%https://community.sw.siemens.com/s/article/what-is-a-power-spectral-density-psd

