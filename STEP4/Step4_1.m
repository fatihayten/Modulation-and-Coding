clear all
clc
close all

F_cut = 1e6;    
F_carrier = 2e9;

% CFO
ppm = 0;
CFO = ppm*F_carrier*1e-6;
random_phase = 0;

F_sampling = 8e6; %sampling frequency
T_sampling = 1/F_sampling;
F_symbol = 2*F_cut; %symbol frequency
T_symbol = 1/F_symbol;

% FILTER
upsample_rate = F_sampling/F_symbol;
beta=0.9;
span = 5; 
taps = upsample_rate*span+1;

% Simulation Parameters
modulation_type = 'qam';
Nbps = 2;
snr = 0:2:16; % We are free to simulate at any SNR
number_of_bits = 2000;
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);

% Frame and Frequency Acquisition Parameters
N_values = [20 40 80];
% Since 10,20 and 40 are number of symbols in Page 58-Slides,
%they should be multiplied by Nbps=2.
K_values = [1 8 16];
iteration = 100; % Std deviation is calculated after many iterations. 
time_error = zeros(iteration,length(snr),length(K_values),length(N_values));
freq_error = zeros(iteration,length(snr),length(K_values),length(N_values));

for k = 1:length(K_values)
for n = 1:length(N_values) 
K = K_values(k);
N = N_values(n);
for iter = 1:iteration

Npilot = N;
pilot_bits = ones(Npilot,1);
pilot_position = 200;

pilot_symbols = mapping(pilot_bits,Nbps,modulation_type);
symbol_tx = [transmitted_symbols(1:pilot_position-1);pilot_symbols;transmitted_symbols(pilot_position:end)];

upsampled_symbols = upsample(symbol_tx,upsample_rate);

% Applying first filter
[filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);
    
E_symbol_av = sum(abs(symbol_tx.*symbol_tx))/length(symbol_tx);
Eb = E_symbol_av/Nbps;


for m = 1:length(CFO)
    cfo = CFO(m);
    for j = 1:length(snr)

        %ADD NOISE
        noisy_symbols = add_noise(snr(j),modulation_type,Eb,filtered_symbols);

%         %ADD CFO
%         t1 = (((0:length(noisy_symbols)-1))*T_sampling)';
%         cfo_exponential = exp(1j*(2*pi*(cfo*t1+random_phase)));
%         noisy_symbols = noisy_symbols.*cfo_exponential;

         %  Applying second filter
        filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

        % Eliminating filtered data
        filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);

        % Downsampling
        filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

        % Frame and frequency acquisition
        [est_n,est_cfo] = cfoEstimate(filtered_noisy_symbols,pilot_symbols,T_symbol,K);
        
%         % CFO Compensation
%         t2 = (((0:length(filtered_noisy_symbols)-1))*T_sampling)';
%         cfo_exponential_compensation = exp(-1j*(2*pi*(cfo*t2+random_phase)));
%         compensated_symbols = filtered_noisy_symbols.*cfo_exponential_compensation;

        time_error(iter,j,k,n) = est_n - pilot_position;
        freq_error(iter,j,k,n) = cfo - est_cfo;

    end
end

end
end
end
%% Plot results
% load CFO_K_N.mat
time_error_mean = std(time_error);  
freq_error_mean = std(freq_error*1e6/F_carrier);

figure
plot(snr,time_error_mean(1,:,2,1),'-r');hold on;
plot(snr,time_error_mean(1,:,2,2),'-g');
plot(snr,time_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdev [samples]');
legend('N = 10, K = 8','N = 20, K = 8','N = 40, K = 8');
title('Pilot ToA error as a function of pilot length')
grid on;
ylim([0 200])

figure
plot(snr,time_error_mean(1,:,1,1),'-r');hold on;
plot(snr,time_error_mean(1,:,2,1),'-g');
plot(snr,time_error_mean(1,:,3,1),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdv [samples]');
legend('N = 20, K = 1','N = 20, K = 8','N = 20, K = 16');
title('Pilot ToA error as a function of K')
grid on;
ylim([0 200])

figure
plot(snr,freq_error_mean(1,:,2,1),'-r');hold on;
plot(snr,freq_error_mean(1,:,2,2),'-g');
plot(snr,freq_error_mean(1,:,2,3),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 10, K = 8','N = 20, K = 8','N = 40, K = 8');
title('CFO error as a function of pilot length')
grid on;

figure
plot(snr,freq_error_mean(1,:,1,2),'-r');hold on;
plot(snr,freq_error_mean(1,:,2,2),'-g');
plot(snr,freq_error_mean(1,:,3,2),'-b');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('N = 20, K = 1','N = 20, K = 8','N = 20, K = 16');
title('CFO error as a function of K')
grid on;