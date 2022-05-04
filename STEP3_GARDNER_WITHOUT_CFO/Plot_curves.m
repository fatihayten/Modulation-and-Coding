clc
clear all
close all

snr=10;
time_shift_vector = [1];
x = zeros(length(time_shift_vector),3); 
y = size(x); 


% 4-QAM
for i=1:y(1)
    for j=1:y(2)
        [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr,'qam',2,time_shift_vector(i));
    end
end
ber_qam_4 = (mean(x,2)/480000).';

% % 16-QAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr,'qam',2);
%     end
% end
% ber_qam_16 = (mean(x,2)/480000).';
% 
% % 64-QAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr,'qam',2);
%     end
% end
% ber_qam_64 = (mean(x,2)/480000).';
% 
% % 256-QAM
% for i=1:y(1)
%     for j=1:y(2)
%         [x(i,j),filter_taps_t,filter_taps_f,upsampled_symbols,filtered_symbols] = ber_calculator(snr,'qam',2);
%     end
% end
% ber_qam_256 = (mean(x,2)/480000).';

%% Plotting QAM BER Curves

% figure;
% % Experimental results
% semilogy(snr,ber_qam_4,'-');
% hold on
% semilogy(snr,ber_qam_16,'-');
% semilogy(snr,ber_qam_64,'-');
% semilogy(snr,ber_qam_256,'-');
% % Indicating the BER range
% ylim([1e-5 1]);
% xlabel('E_B/N_0 [dB]');
% ylabel('BER');
% title(['Comparison of Experimental and Theoretical Results of Bit Error',...
%     ' Rate vs E_B/N_0 (dB) for Different Modulation Types']);
% grid on;
% legend('4-QAM Experimental','16-QAM Experimentaml','64-QAM Experimental',...
%     '256-QAM Experimental','4-QAM Theoretical','16-QAM Theoretical',...
%     '64-QAM Theoretical','256-QAM Theoretical')
