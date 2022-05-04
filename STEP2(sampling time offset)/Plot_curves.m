clc
clear all
close all

snr=(-5:1:15);
x = zeros(length(snr),1); % Number indicates the simulation repeat number, 1 is ok, 3 is ideal.
y = size(x); 
shift = linspace(0,0.1,5); % 1,2,3,4 and 5 sampling time shifts will be invstigated.
%% QAM Calculations
% 1 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4,1);
    end
end
ber_qam_1 = (mean(x,2)/480000).';

% 2 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4,2);
    end
end
ber_qam_2 = (mean(x,2)/480000).';

% 3 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4,3);
    end
end
ber_qam_3 = (mean(x,2)/480000).';

% 4 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4,4);
    end
end
ber_qam_4 = (mean(x,2)/480000).';

% 5 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr(i),'qam',4,5);
    end
end
ber_qam_5 = (mean(x,2)/480000).';

%% Plotting QAM BER Curves

figure;
% Experimental results
semilogy(snr,ber_qam_1,'--');
hold on
semilogy(snr,ber_qam_2,'-o');
semilogy(snr,ber_qam_3,'-+');
semilogy(snr,ber_qam_4,'-*');
semilogy(snr,ber_qam_5,'-*');
ylim([1e-3 1]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER degradation for increasing values of the sample time shift');
grid on;
legend('no time offset','t_0=0.025 T','t_0=0.05 T',...
    't_0=0.075 T','t_0=0.1 T')