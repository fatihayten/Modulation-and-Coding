function [error_symbols,filter_tapss,filter_tapss_freq,upsampled_symbols,filtered_symbols,error] = ber_calculator(snr,modulation_type,Nbps,shift)

number_of_bits = 480000;
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);


cfo = 0;
K = 0.07;
fc = 2e+9; % carrier frquency
f_cut = 1e6; % cutoff frquency
F_sampling = 80e6; %sampling frequency, it has been increase significantly to 
% be able to simulate time shift
T_sampling = 1/F_sampling;
F_symbol =2*f_cut; %symbol frequency
T_symbol = 1/F_symbol;
upsample_rate = F_sampling/F_symbol;

upsampled_symbols = upsample(transmitted_symbols,upsample_rate);

beta=0.3;
span = 10;
taps = upsample_rate*span+1;

% Applying first filter
[filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);


E_symbol_av = sum(abs(transmitted_symbols.*transmitted_symbols))/length(transmitted_symbols);

%Energy Per Bit, To find it, I divide symbol energy into bits per symbol
Eb = E_symbol_av/Nbps;

% Adding noise to symbols
noisy_symbols = add_noise(snr,modulation_type,Eb,filtered_symbols);

% %Adding CFO
% t1 = ((0:length(noisy_symbols)-1))*(1/(F_sampling));
%     
% for i = 1:length(t1)
%     noisy_symbols(i) = noisy_symbols(i)*exp(1j*(2*pi*cfo*t1(i)));
% end

%  Applying second filter
filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

% Eliminating filtered data
filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);

%CIRCLING
% With this function i am circling the matrix, i.e., [1 2 3 4] becomes [2
% 3 4 1] if it is circled 1 times.
filtered_noisy_symbols = circshift(filtered_noisy_symbols,shift);
% In our case Tsymbol = 1/2e6 and Tsampling = 1/80e6. To simulate 0.1Tsymbol time
% shift, the required sample shift = 0.1*(1/2e6) / (1/80e6) = 4.
% Likewise,
% 1 sample shift >= 0.025 Tsymbol shift
% 2 sample shift >= 0.05 Tsymbol shift
% 3 sample shift >= 0.075 Tsymbol shift
% 4 sample shift >= 0.1 Tsymbol shift

% %Compensate CFO
% t2 = ((0:length(filtered_noisy_symbols)-1))*(1/(F_sampling));
%     
% for i = 1:length(t2)
%     filtered_noisy_symbols(i) = filtered_noisy_symbols(i)*exp(1j*(2*pi*cfo*t2(i)));
% end

    % GARDNER
   
    M=upsample_rate;
    symbol_rx_upsampled = filtered_noisy_symbols;
    L = length(symbol_rx_upsampled);
    L = L-mod(L,M);
            
    error = zeros(L/M,1);
    corr = zeros(L/M,1);
    prevY = symbol_rx_upsampled(1);
            
    for i = 1:(L/M)-1
        a = ((i-1)*M:M*i-1);
        b = symbol_rx_upsampled(1+(i-1)*M:i*M);
        c = M/2+(i-1)*M-error(i);
        c2 = i*M-error(i);
        
        Y_mid = interp1(a,b,c,'pchip');
        Y = interp1(a,b,c2,'pchip');
        
        corr(i) = (2*K)*real(Y_mid*(conj(Y) - conj(prevY)));
        error(i+1) = error(i) + corr(i);
        prevY = Y;
    end
    
    % Error Correction
    %timeshift = round((shift-error(end)));
    filtered_noisy_symbols = circshift(filtered_noisy_symbols,round(error(end)));

% Downsampling
filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

% Demapping
demodulated_bits = demapping(filtered_noisy_symbols,Nbps,modulation_type);


error_symbols = sum(bits ~= demodulated_bits);


filter_tapss = filter_taps_in_time;
filter_tapss_freq=filter_taps_in_freq;
end


