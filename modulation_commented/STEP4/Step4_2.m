clear all
clc
close all

%% Simulation Parameters
F_cut = 1e6;    
F_carrier = 2e9;
F_sampling = 400e6; %sampling frequency
T_sampling = 1/F_sampling;
F_symbol = 2*F_cut; %symbol frequency
T_symbol = 1/F_symbol;
modulation_type = 'qam';
Nbps = 2;
snr = 0:2:16; % We are free to simulate at any SNR
number_of_bits = 400;
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);

shift = 4; % Sampling offset to simulate 0.1T_symbol

%% CFO
ppm = [0 10];
CFO = ppm*F_carrier*1e-6;

%% Filter Parameters
upsample_rate = F_sampling/F_symbol;
beta=0.3;
span = 1; 
taps = upsample_rate*span+1;

%% Frame and Frequency Acquisition Parameters
N_values = 80; %corresponds to 40 because 
K_values = 8;
iteration = 200; % Std deviation is calculated after many iterations. 
% Decrease it if it takes long time


time_error = zeros(iteration,length(CFO),length(snr));
freq_error = zeros(iteration,length(CFO),length(snr));

 
K = K_values;
N = N_values;


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

for iter = 1:iteration
    for m = 1:length(CFO)
    cfo = CFO(m);
        for j = 1:length(snr)
        
        %ADD NOISE
        noisy_symbols = add_noise(snr(j),modulation_type,Eb,filtered_symbols);

        %ADD CFO
        t1 = (((0:length(noisy_symbols)-1))*T_sampling)';
        cfo_exponential = exp(1j*(2*pi*(cfo*t1)));
        noisy_symbols = noisy_symbols.*cfo_exponential;

         %  Applying second filter
        filtered_noisy_symbols = conv(noisy_symbols,filter_taps_in_time);

        % Eliminating filtered data
        filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);

        if cfo == 10
            %Sampling Time Offset
            filtered_noisy_symbols = circshift(filtered_noisy_symbols,shift);
    
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
    
    
            %Sampling Time Offset Error Correction
            filtered_noisy_symbols = circshift(filtered_noisy_symbols,round(error(end)));
        
        end

        % Downsampling
        filtered_noisy_symbols = downsample(filtered_noisy_symbols,upsample_rate);

        % Frame and frequency acquisition
        [est_n,est_cfo] = acquisition(filtered_noisy_symbols,pilot_symbols,T_symbol,K);
        
        % CFO Compensation
        t2 = (((0:length(filtered_noisy_symbols)-1))*T_sampling)';
        cfo_exponential_compensation = exp(-1j*(2*pi*(est_cfo*t2)));
        compensated_symbols = filtered_noisy_symbols.*cfo_exponential_compensation;

        time_error(iter,m,j) = est_n - pilot_position;
        freq_error(iter,m,j) = cfo - est_cfo;

        end
    end

end

%% Plot results
% load CFO_K_N.mat
time_error_mean = std(time_error);  
freq_error_mean = std(freq_error*1e6/F_carrier);

figure
plot(snr,reshape(time_error_mean(1,1,:),[1,length(snr)]),'-r');hold on;
plot(snr,reshape(time_error_mean(1,2,:),[1,length(snr)]),'-g');
xlabel('E_B/N_0 [dB]');
ylabel('Time error stdev [samples]');
legend('no CFO, \epsilon=0','10 ppm CFO, \epsilon=0.02');
title('Pilot ToA error if CFO and time error,N = 40, K = 8')
grid on;
ylim([0 2])


figure
plot(snr,reshape(freq_error_mean(1,1,:),[1,length(snr)]),'-r');hold on;
plot(snr,reshape(freq_error_mean(1,2,:),[1,length(snr)]),'-g');
xlabel('E_B/N_0 [dB]');
ylabel('Frequency error stdv [ppm]');
legend('no CFO, \epsilon=0','10 ppm CFO, \epsilon=0.02');
title('CFO error if CFO and time error,N = 40, K = 8')
grid on;
ylim([0 2])
