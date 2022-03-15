function [noisy_symbols] = add_noise(snr,modulation_type,Eb,filtered_symbols)
% If the modulation type is PAM, Gaussian noise will be added in 1D fashion
% If the modulation type is QAM, Gaussian noise will be added in 2D fashion

% IN PROJECT, WE ONLY CONSIDER QAM MODULATION
% if modulation_type == "pam"
%     for j = 1:length(snr)
%     N0 = Eb/(10^(snr/10));
%     NoisePower = N0;
%     noise = sqrt(NoisePower/2)*(randn(length(filtered_symbols),1));
%     end
% 
% end

% Since we know the Eb, we should find N0.
% After finding N0, we have to divide this power into two parts;
% real and imaginary axises. Because QAM symbols exist in complex plane.
% Power in real axis will be N0/2; power in imag axis also will be N0/2.

% randn function produces zero mean unit(1) variance gaussian distributed
% numbers. Variance for gaussian variables is very crucial because it determines
% the power. So, in the beginning (when we write randn(...)), power = variance = 1.
% When we multiply a random variable with constant "a", its variance becomes
% a*a.
% To reach N0/2 power (variance) in each axis, we are multiplying the random variables by 
% sqrt (N0/2).

% We are producing noise for different SNR values, using a loop.
if modulation_type == "qam"
    for k = 1:length(snr)
    N0 = Eb/(10^(snr/10));
    NoisePower = N0;
    noise = sqrt(NoisePower/2)*(randn(length(filtered_symbols),1)+1i*randn(length(filtered_symbols),1));
    end
end

% Noise is added to symbols
noisy_symbols = noise + filtered_symbols; 
end