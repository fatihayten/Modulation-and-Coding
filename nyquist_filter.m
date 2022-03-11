function [result] = nyquist_filter(beta,span,sps,F_samp,F_symb)

% To create this filter, we get inspired from the "rcosdesign" function.
% Ideal raised cosine filters have an infinite number of taps (i.e., has 
% infinite impulse response). 
% Therefore, practical raised cosine filters are windowed. 
% The window length is controlled using the span and sps arguments.
% Span is number of symbols; and sps is the number of samples per symbol 
% (oversampling factor
% For example, if span is specified as 4, filter length will be equal to 4
% symbol duration.
% Also, since the symbols are already upsampled, the number of taps will be
% equal to span*sps+1 (filter order = span*sps).
% The order of the filter must be even.

% Beta is roll-off factor of the filter. It determines the excess bandwidth of the
% filter. To implement a meaningful filter, beta should be in range [0,1].

if beta<0 | beta >1
    error("Beta(roll-off) parameter of Nyquist Filter should be in range [0,1]");
end
T = (1/F_samp)*sps; % symbol duration, = F_sybm
taps = span*sps+1;

stepo = (1 / taps) * F_samp;
highf = stepo * (taps - 1)/2;
f = linspace(-highf, highf, taps);
lcorner = (1 - beta) / (2 * T);
rcorner = (1 + beta) / (2 * T);
result = ones(1, taps) .* (sqrt(T / 2 * (1 + cos(pi * T / beta * (abs(f) - lcorner)))));
result(abs(f) < lcorner) = sqrt(T);
result(abs(f) > rcorner) = 0;
figure
plot(linspace(-F_samp/2,F_samp/2,taps),result)
title("Fourier Transform of Nyquist Filter")
xlabel("Frquency (Hz)")
result_t = fftshift(ifft(ifftshift(result)));


energy = sqrt(sum(result_t.*result_t));
result = result_t/energy; % To implement unit energy filter


end