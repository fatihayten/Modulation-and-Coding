function [error_symbols,filter_tapss,filter_tapss_freq,upsampled_symbols,filtered_symbols] = ber_calculator(snr,modulation_type,Nbps)
% This function is created for calculating Bit Error Rate.

number_of_bits = 480000; % Ideal : 480000
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);

% You can observe the "SYMBOLS" using scatterplot function.
%scatterplot(transmitted_symbols,10) 

F_sampling = 8e6; %sampling frequency
F_symbol = 2e6; %symbol frequency
upsample_rate = F_sampling/F_symbol; %Should be integer, because L argument of
% upsample function should be integer-valued.
% https://nl.mathworks.com/help/control/ref/lti.upsample.html?s_tid=doc_ta


upsampled_symbols = upsample(transmitted_symbols,upsample_rate);

beta=0.3;
span = 10; %I prefer 10, works perfectly, you can try other span values
taps = upsample_rate*span+1; % number of taps (number of samples in time domain)

% Applying first filter
[filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);

%  Applying first filter
filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);

% Energy per Symbol, Reference: "Digital Modulation" Slideset, Page 14
% Basically, I sum all of energy of symbols and divide it to number of
% symbols in order to find the "average" symbol energy.
E_symbol_av = sum(abs(transmitted_symbols.*transmitted_symbols))/length(transmitted_symbols);

%Energy Per Bit, To find it, I divide symbol energy into bits per symbol
Eb = E_symbol_av/Nbps;

% Adding noise to symbols
noisy_symbols = add_noise(snr,modulation_type,Eb,filtered_symbols);

%  Applying second filter
filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

% Eliminating filtered data
filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);
%I couldn understand completely, but i guess, when we do convolution
%operation, length of the result  is M+N-1 (from DSP). In real life
% conv operation is done by FTT and IFFT which gives the result that has same
% length as input. (check conv(x,y,"same")). Therefore, we are eliminating
% extra values and only consider the data in the middle of the result.

% Downsampling
filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

% Demapping
demodulated_bits = demapping(filtered_noisy_symbols,Nbps,modulation_type);

% Finding The Probability of Error in Bits (by checking the matrices
% elementwise)
error_symbols = sum(bits ~= demodulated_bits);

% Following 2 lines are not important, just about input-output issues of
% the function
filter_tapss = filter_taps_in_time;
filter_tapss_freq=filter_taps_in_freq;
end
%Power Spectral Density, A Good Resource
%https://community.sw.siemens.com/s/article/what-is-a-power-spectral-density-psd

