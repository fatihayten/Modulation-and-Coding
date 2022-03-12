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

% Beta is roll-off factor of the filter. It determines the excess bandwidth...
% of the filter. To implement a meaningful filter, beta should be in range...
% [0,1].

if beta<0 || beta >1
    error("Beta(roll-off) parameter of Nyquist Filter should be in " + ...
        "range [0,1]");
end
T = (1/F_samp)*sps; % symbol duration, = F_sybm
taps = span*sps+1;

% Using the frequency component "F_sample/2"
% is not wise, therefore, we leave a margin from positive and negative
% parts of the frequency axis. Also, sthis margin should be dependent on
% taps and sampling frequency.
freq_margin = (1/taps)*F_samp; 

% Maximum frequency component of the filter
highf = freq_margin*(taps-1)/2; 

% Dividing the frequency axis into pieces
f = linspace(-highf, highf, taps); 

% Deciding the frequency limits according to frequency response formula of
% Nyquist Filters ("Digital Communications, Fourth Edition", J. G. Proakis,
% 4th Edition, Page 546)
limit1 = (1-beta)/(2*T);
limit2 = (1+beta)/(2*T);

% Since we will divide the filter into two parts, we should use the formula
% after applying square root operation.
filter_in_freq_domain = ones(1, taps).*...
    (sqrt(T/2*(1+cos(pi*T/beta*(abs(f)-limit1)))));
filter_in_freq_domain(abs(f)<limit1)=sqrt(T);
filter_in_freq_domain(abs(f)>limit2) = 0;

% figure;
% plot(linspace(-F_samp/2,F_samp/2,taps),filter_in_freq_domain)
% title("Fourier Transform of Nyquist Filter")
% xlabel("Frquency (Hz)")

% Transforming the filter into time domain
filter_in_time_domain = fftshift(ifft(ifftshift(filter_in_freq_domain)));

% To implement unit energy filter
result = filter_in_time_domain / sqrt(sum(filter_in_time_domain.^2)); 


end