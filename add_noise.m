function [noisy_symbols] = add_noise(snr_range,modulation_type,Eb,filtered_symbols,F_sampling)
% If the modulation type is PAM, Gaussian noise will be added in 1D fashion
% If the modulation type is QAM, Gaussian noise will be added in 2D fashion

if modulation_type == "pam"
    for j = 1:length(snr_range)
    N0 = Eb/10.^(snr_range(j)/10);
    NoisePower = N0;
    noise = sqrt(NoisePower)*(randn(length(filtered_symbols),1));
    end
end
if modulation_type == "qam"
    for k = 1:length(snr_range)
    N0 = Eb/10.^(snr_range(k)/10);
    NoisePower = N0;
    noise = sqrt(NoisePower/2)*(randn(length(filtered_symbols),1)+1i*randn(length(filtered_symbols),1));
    end
end
noisy_symbols = noise + filtered_symbols;
end