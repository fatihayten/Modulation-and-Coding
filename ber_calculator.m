function [error_symbols,filter_tapss,filter_tapss_freq,upsampled_symbols,filtered_symbols] = ber_calculator(snr,modulation_type,Nbps)
 
number_of_bits = 4800; % IDeal : 480000
bits = random_bit_generator(number_of_bits);
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
span = 40; %IDEAL:10
taps = upsample_rate*span+1;
[filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
% figure;
% plot(a)
% hold on

%checking the difference between self-written filter and rcosdesign filter
% b = rcosdesign(beta,span,upsample_rate,"normal");
% plot(b)
% legend("Self-written filter","rcosdesing filter")
% difference_of_2_filters = sum(abs(a-b));
% title("Nyquist Filter in Time Domain")

%illustration of cancellation of the inter-symbol interference
%  figure;
%  hold on
%  plot(conv(filter_taps,filter_taps))
%  plot(circshift(conv(filter_taps,filter_taps),upsample_rate))
%  plot(circshift(conv(filter_taps,filter_taps),2*upsample_rate))


% Passing the symbols through the filter1
filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);


% Average Symbol Energy
E_symbol_av = sum(abs(transmitted_symbols))/length(transmitted_symbols);

%Energy Per Bit
Eb = E_symbol_av/Nbps;



% Adding noise to symbols
noisy_symbols = add_noise(snr,modulation_type,Eb,filtered_symbols);

% Passing the symbols through the filter1
filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

% Eliminating filtered data
filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);

% Downsampling
filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

% Demapping
demodulated_bits = demapping(filtered_noisy_symbols,Nbps,modulation_type);

% Finding The Probability of Error in Bits
error_symbols = sum(bits ~= demodulated_bits);
filter_tapss = filter_taps_in_time;
filter_tapss_freq=filter_taps_in_freq;
end
%https://community.sw.siemens.com/s/article/what-is-a-power-spectral-density-psd

