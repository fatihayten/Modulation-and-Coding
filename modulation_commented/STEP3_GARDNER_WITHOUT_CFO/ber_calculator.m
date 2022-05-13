clear all
clc
close all

K_values = [0.01 0.03 0.05 0.2]; % Gardner Parameters

%% Simulation Parameters
modulation_type = 'qam';
Nbps = 2;
snr = 15; % We are free to simulate at any SNR
number_of_bits = 4800;
bits = random_bit_generator(number_of_bits);
transmitted_symbols = mapping(bits,Nbps,modulation_type);
F_cut = 1e6;    
F_sampling = 80e6; %sampling frequency
F_symbol = 2*F_cut; %symbol frequency
T_symbol = 1/F_symbol;
upsample_rate = F_sampling/F_symbol; %Should be integer, because L argument of 
upsampled_symbols = upsample(transmitted_symbols,upsample_rate);
 
%% Filter Parameters
beta=0.3;
span = 5; 
taps = upsample_rate*span+1;

%% Time Shift
shift = [1 2 3 4]; % Shift should be between 1 and 40 since upsample rate is 40.

%% Simulation
iteration = 10; % To get average error, we will simulate multiple times.
for x = 1:length(shift)
for m = 1:length(K_values)
K = K_values(m);
for iter = 1:iteration

    % Applying first filter
    [filter_taps_in_time,filter_taps_in_freq] = nyquist_filter(beta,span,upsample_rate,F_sampling,F_symbol);
    
    %  Applying first filter
    filtered_symbols = conv(upsampled_symbols,filter_taps_in_time);
    
    
    E_symbol_av = sum(abs(transmitted_symbols.*transmitted_symbols))/length(transmitted_symbols);
    Eb = E_symbol_av/Nbps;
    
    % Adding noise to symbols (NO NOISE IN THIS PART)
    %noisy_symbols = add_noise(snr,modulation_type,Eb,filtered_symbols);
    
    %  Applying second filter
    filtered_noisy_symbols = conv(filtered_symbols,filter_taps_in_time);
   
    % Eliminating filtered data
    filtered_noisy_symbols = filtered_noisy_symbols(taps:end-taps+1);
    
    % Sampling Time Error
    % this time circshift method didnt work, so I utilized indexing to simulate time shift offset
    filtered_noisy_symbols = filtered_noisy_symbols(1+shift(x):end);

    % GARDNER Algorithm
    M=upsample_rate;
    L = length(filtered_noisy_symbols);
    L = L-mod(L,M);
            
    error = zeros(L/M,1);
    corr = zeros(L/M,1);
    initialization = filtered_noisy_symbols(1);
            
    for i = 1:(L/M)-1
        a = ((i-1)*M:M*i-1);
        b = filtered_noisy_symbols(1+(i-1)*M:i*M);
        c = M/2+(i-1)*M-error(i);
        c2 = i*M-error(i);
        
        Y_mid = interp1(a,b,c,'pchip');
        Y = interp1(a,b,c2,'pchip');
        
        corr(i) = (2*K)*real(Y_mid*(conj(Y) - conj(initialization)));
        error(i+1) = error(i) + corr(i);
        initialization = Y;
    end
    
    % Since our aim is to 'find' and 'track' the error, we want to keep the difference
    % between real error (variable 'shift') and our estimation (variable 'error')
    shift_vector = ones(size(error))*shift(x);
    time_error(iter,:,m,x) = (shift_vector-error)'*T_symbol;

       
    % No need to down sampling and demapping
end
end


% %% By using 'smooth' function, final waveforms will be put in a better shape.
% 
% for m = 1:length(K_values)
% for iter = 1:iteration
%     time_error(iter,:,m,x) = smooth(time_error(iter,:,m,x));
% 
% end
% end

%% Plotting
time_error_mean = mean(time_error);
time_error_stdv = std(time_error);
mean1 = time_error_mean(1,:,1,x);
mean2 = time_error_mean(1,:,2,x);
mean3 = time_error_mean(1,:,3,x);
mean4 = time_error_mean(1,:,4,x);
stdv1 = time_error_stdv(1,:,1,x);
stdv2 = time_error_stdv(1,:,2,x);
stdv3 = time_error_stdv(1,:,3,x);
stdv4 = time_error_stdv(1,:,4,x);

subplot(2,2,x);
plot(mean1,'r','LineWidth',1);hold on;
plot(mean2,'g','LineWidth',1)
plot(mean3,'b','LineWidth',1)
plot(mean4,'c','LineWidth',1)
plot(mean1+stdv1,'--r')
plot(mean1-stdv1,'--r')
plot(mean2+stdv2,'--g')
plot(mean2-stdv2,'--g')
plot(mean3+stdv3,'--b')
plot(mean3-stdv3,'--b')
plot(mean4+stdv4,'--c')
plot(mean4-stdv4,'--c')
ylim([-1e-6 3e-6])
xlabel('Symbols');
ylabel('Time error (mean \pm stdv)');
legend('K = 0.01','K = 0.03','K = 0.05','K = 0.2');
a = 0.025*x;
title(sprintf('Convergence of the Gardner algorithm When Sample Shift Error = %.3fT_s_y_m_b_o_l',x*0.025))
grid on;
end