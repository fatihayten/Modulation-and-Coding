clc
clear all
close all

snr=(-5:1:15); % SNR range
x = zeros(length(snr),2); % Number indicates the simulation repeat number, 1 is ok, 3 is ideal.
y = size(x); 

%% CFO = 0 
for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,0);
    end
end
ber_qam_1 = (mean(x,2)/480000).';

%% CFO = 2
for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,2);
    end
end
ber_qam_2 = (mean(x,2)/480000).';

%% CFO = 10 
for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,10);
    end
end
ber_qam_3 = (mean(x,2)/480000).';


%% Plotting BER Curves

figure;
% Experimental results
ber_qam_4_theoretical = berawgn(snr,'qam',16); % 4-QAM BER Curve
% semilogy(snr,ber_qam_4_theoretical,'--');


semilogy(snr,ber_qam_1,'-*');
hold on
semilogy(snr,ber_qam_2,'-o');
semilogy(snr,ber_qam_3,'-+');
 ylim([1e-5 1]);
 xlim([-5 15]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('Impact of CFO on BER (only ISI)');
grid on;
legend('no CFO','CFO = 2 ppm',...
    'CFO = 10 ppm')