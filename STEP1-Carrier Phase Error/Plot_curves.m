clc
clear all
close all

snr=(-5:1:15);
x = zeros(length(snr),2); % Number indicates the simulation repeat number, 1 is ok, 3 is ideal.
y = size(x); 
%% QAM Calculations
% 1 sample time shift

for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,0);
    end
end
ber_qam_0 = (mean(x,2)/1e6).';



for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,1/120);
    end
end
ber_qam_3 = (mean(x,2)/1e6).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,1/72);
    end
end
ber_qam_4 = (mean(x,2)/1e6).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,1/36);
    end
end
ber_qam_5 = (mean(x,2)/1e6).';

for i=1:y(1)
    for j=1:y(2)
        [x(i,j)] = ber_calculator(snr(i),'qam',4,1/12);
    end
end
ber_qam_6 = (mean(x,2)/1e6).';


%% Plotting QAM BER Curves

figure;
% Experimental results
ber_qam_4_theoretical = berawgn(snr,'qam',16); % 4-QAM BER Curve
% semilogy(snr,ber_qam_4_theoretical,'--');

semilogy(snr,ber_qam_0,'-*');
hold on
semilogy(snr,ber_qam_3,'-*');
semilogy(snr,ber_qam_4,'-*');
semilogy(snr,ber_qam_5,'-*');
semilogy(snr,ber_qam_6,'-*');

%semilogy(snr,ber_qam_4,'');
 ylim([1e-5 1]);
 xlim([-5 15]);
xlabel('E_B/N_0 [dB]');
ylabel('BER');
title('Impact of Carrier Phase Error on BER');
grid on;
legend("no phase error ", "phase error = 3" + char(176) , "phase error = 5" + char(176) , "phase error = 10" + char(176) , "phase error = 30" + char(176),'Location','best')


