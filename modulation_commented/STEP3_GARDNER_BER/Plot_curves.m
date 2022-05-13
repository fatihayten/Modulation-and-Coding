clc
clear all
close all

snr=(-5:1:15); % SNR Range
x = zeros(length(snr),2); % Number indicates the simulation repeat number, 1 is ok, 3 is ideal.
y = size(x); 
%% QAM Calculations

%% 1 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols,error1] = ber_calculator(snr(i),'qam',4,1);
    end
end
ber_qam_1 = (mean(x,2)/480000).';

%% 2 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols,error2] = ber_calculator(snr(i),'qam',4,2);
    end
end
ber_qam_2 = (mean(x,2)/480000).';

%% 3 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols,error3] = ber_calculator(snr(i),'qam',4,3);
    end
end
ber_qam_3 = (mean(x,2)/480000).';

%% 4 sample time shift
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols,error4] = ber_calculator(snr(i),'qam',4,4);
    end
end
ber_qam_4 = (mean(x,2)/480000).';


%% Plotting QAM BER Curves

figure;
% Experimental results
ber_qam_4_theoretical = berawgn(snr,'qam',16); % 4-QAM BER Curve
semilogy(snr,ber_qam_4_theoretical,'--');
hold on
semilogy(snr,ber_qam_1,'gs');
semilogy(snr,ber_qam_2,'-o');
semilogy(snr,ber_qam_3,'-+');
semilogy(snr,ber_qam_4,'-*');
ylim([1e-5 1]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('BER degradation for increasing values of the sample time shift after Gardner Algorithm is Employed');
grid on;
legend('Theoretical','t_0=0.025 T','t_0=0.05 T',...
    't_0=0.075 T','t_0=0.1 T')