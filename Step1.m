number_of_bits = 100;
bits = random_bit_generator(number_of_bits);
Nbps = 2;
modulation = "qam";
transmitted_symbols = mapping(bits,Nbps,modulation);
%scatterplot(transmitted_symbols,10)
%plot(transmitted_symbols,'g*')
F_sampling = 8e6; %sampling frequency
F_symbol = 2e6; %symbol frequency
upsample_rate = F_sampling/F_symbol; %Should be integer, because L argument of...
% upsample function should be integer-valued.
% https://nl.mathworks.com/help/control/ref/lti.upsample.html?s_tid=doc_ta

upsampled_symbols = upsample(transmitted_symbols,upsample_rate);

beta=0.3;
span = 25;
a = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol)
figure
plot(a)
title("Nyquist Filter in Time Domain")