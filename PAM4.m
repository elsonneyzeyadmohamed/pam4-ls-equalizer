clc ;
clear ;
close all;
Fs = 40000;
Ts = 1/Fs;
%% rect pulse 
rect_length = 20 ;

Tb = rect_length * Ts ; % symbol duration

main_lobe_freq = 1/Tb;

rect = ones(1, rect_length) ; % here we start to declare the bit duration which is 6 ;
random_bits_N = 20000; % bits/2 PAM4 symbols 
bits = randi([0, 1], 1, random_bits_N);

% 2 bits -> 1 PAM-4 symbol
bits_2 = reshape(bits, 2, []).';
N_symbols = length(bits_2);
symbol_values = bits_2(:,1)*2 + bits_2(:,2);
symbol_values = symbol_values(:).';


%% signal generation 
% PAM-4 mapping
signal_values = -2 + (4/3)*symbol_values;

signal = upsample(signal_values , rect_length); % insert zeros to leave space for convolution to take place for the hold
PAM_conv = conv(rect, signal );
% To truncate the final values that were generated from convolution
PAM = PAM_conv(1: end - rect_length + 1) ;
N = length(PAM);
t_observ1 = N*Ts;

time1 = 0 : Ts : t_observ1 - Ts ;

%% freq domain
% for rect
N_fft = N*100;
PAM_fd1 =  fftshift(fft(PAM, N_fft)) *Ts;
f1 = 1/(N_fft*Ts);
freq_domain = (-N_fft/2:1:N_fft/2 - 1) * f1; 
%% single pole channel

fc_sing_pole = main_lobe_freq/4 ;

s = freq_domain * 1j;

tran_fun = 1./(1 + s/fc_sing_pole);
single_pole_channel_t = ifft(ifftshift(tran_fun))/Ts;
%% Applying single pole channel

single_pole_out = tran_fun .* PAM_fd1;
PAM_sing_pole = ifft(ifftshift(single_pole_out))/Ts;
%PAM_sing_pole = align_signals(PAM_sing_pole_full, PAM);
PAM_sing_pole = PAM_sing_pole(1 : N);
% frequency :
single_pole_fd = fftshift(fft(PAM_sing_pole, N_fft)) *Ts;
%% Ideal low-pass filter
filter_BW = main_lobe_freq ;
low_pass_filter_fd = (freq_domain < filter_BW & (freq_domain >= (-filter_BW)) );
low_pass_filter = ifft(ifftshift(low_pass_filter_fd), 'symmetric')/Ts;

%% Apply receiver LPF to clean signal
filtered_signal_fd_clean = single_pole_fd .* low_pass_filter_fd;
filtered_signal_clean = ifft(ifftshift(filtered_signal_fd_clean), 'symmetric')/Ts;
filtered_signal_clean = filtered_signal_clean(1 : N);
%% LS Equalizer

pulse = filtered_signal_clean;    %brickwall_out ;     % Received signal
N = length(pulse);
OS = rect_length;             % Oversampling factor (20)
L = 8;                      % Equalizer length

MSE = zeros(L , 1);   % error calculation 
%% Try all possible target delays
% Desired original signal
    target_bits_values =  pulse;
    
best_MSE = inf;
best_tau = 1;
best_i = 1;
best_w = zeros(L, 1);

for k = 1:OS
    downsampled_target_sig = pulse(k : OS : OS*N_symbols);   % the instance point 
    downsampled_target_sig = downsampled_target_sig(:);  % turn to column matrix 
   
    P = zeros(N_symbols + L - 1, L); % convolution matrix
    % assign values in the matrix
    for col = 1:L
        P(col : col + N_symbols - 1, col) = downsampled_target_sig;   
    end
    P_inv = (P' * P) \ P';

% inner loop to get the right target which sees which is better to care
% about precursors , or post cursors or both 

    for i = 1:L        % L delays only 

    % Target signal for this trial
    Ti = zeros(N_symbols + L - 1, 1);

    
    % Shift target according to the trial delay
    start_idx = i;
    end_idx = i + N_symbols - 1;
    Ti(start_idx:end_idx) = signal_values(:);  % if the matrix depth is finished do not use the rest 

    % solution
    w = P_inv * Ti;

    % Equalizer output to compare
    y = P * w;
    valid_y = y(start_idx:end_idx);
    error = valid_y - signal_values(:);

    % MSE
    MSE(i) = mean((abs(error)).^2);
    if MSE(i) < best_MSE
            best_MSE = MSE(i);
            best_tau = k;
            best_i = i;
            best_w = w;
            best_P = P; % Save the winning matrix for final reconstruction
            best_MSE_array = MSE; %Save the winning curve
    end
    end % end of inner for loop 

end 
%% Find best target delay
fprintf('Best Sampling Phase (tau) = %d\n', best_tau);
fprintf('Best Equalizer Target Index = %d\n', best_i);
fprintf('Absolute Minimum MSE = %.6f\n\n', best_MSE);
%% fixed equalizer for all the coming 
upsampled_equalizer = upsample(best_w,rect_length, best_tau - 1); 
equalizer_fd = fftshift(fft(upsampled_equalizer, N_fft));
equalizer_fd = equalizer_fd(:).';
combined_noise_transfer = low_pass_filter_fd .* equalizer_fd;
% noise enhancement factor
noise_enhancement_factor = mean(abs(combined_noise_transfer).^2); 
%% receiver output with no noise to get isi power 
equalized_output_raw_clean = conv(upsampled_equalizer, filtered_signal_clean);
% exact deterministic alignment 
eq_delay = (best_i - 1) * rect_length + best_tau - 1;
equalized_output_filtered_clean = equalized_output_raw_clean(1 + eq_delay : length(PAM) + eq_delay);
%% Decision sampling
sampling_offset = best_tau; 
clean_sampled_values = equalized_output_filtered_clean(sampling_offset:rect_length:end);
clean_sampled_values = clean_sampled_values(:);
%% Effective signal gain
signal_power = mean(abs(signal_values).^2);
effective_gain = mean(clean_sampled_values .* signal_values(:)) / signal_power; %power difference for isi
%% ISI after complete receiver
residual_error = clean_sampled_values - effective_gain*signal_values(:);
residual_ISI_power = mean(abs(residual_error).^2);
%% Receiver information
fprintf('--- Receiver Analysis ---------\n');
fprintf('Transmitted Signal Power = %.6f\n', signal_power);
fprintf('Effective Equalizer Gain = %.6f\n', effective_gain);
fprintf('Residual ISI Power = %.6f\n', residual_ISI_power);
fprintf('Noise Enhancement Factor = %.6f\n\n', noise_enhancement_factor);

%% Noise addition
EbN0_dB = 0:2:24;                     % The range of Eb/No values to test
BER_simulated = zeros(1, length(EbN0_dB));
BER_theoretical = zeros(1, length(EbN0_dB));


for idx = 1:length(EbN0_dB)

current_EbN0_dB = EbN0_dB(idx);
SNR_linear = 10^(current_EbN0_dB / 10);

% Noise variance for single-sample decision
theoretical_signal_power = mean(abs(signal_values).^2);  % 1 , no ISI
noise_power = (theoretical_signal_power * OS) / (2 * SNR_linear);
noise = sqrt(noise_power) * randn(size(PAM_sing_pole));
noisy_sig =  PAM_sing_pole + noise; 
% frequency domain of noisy signal 
noisy_sig_fd =  fftshift(fft(noisy_sig, N_fft)) *Ts;
% noisy signal filtered 
filtered_signal_fd = noisy_sig_fd .* low_pass_filter_fd;
filtered_signal = ifft(ifftshift(filtered_signal_fd),'symmetric')/Ts;
% Equalizer
equalized_output_raw_filt = conv(upsampled_equalizer,filtered_signal);
equalized_output_filtered = equalized_output_raw_filt(1 + eq_delay : length(PAM) + eq_delay);
equalized_output_fd = fftshift(fft(equalized_output_filtered, N_fft));
if current_EbN0_dB == 0
        saved_very_noisy_equalized = equalized_output_filtered;
end

if current_EbN0_dB == 18
        saved_very_noisy_equalized2 = equalized_output_filtered;
end
%% Comparator
received_symbols = zeros(1, N_symbols);
%sampling_offset = round(rect_length / 2);
sampled_values = equalized_output_filtered(sampling_offset : rect_length : end);

received_symbols(sampled_values >= 4/3) = 3;
received_symbols(sampled_values >= 0 & sampled_values < 4/3) = 2;
received_symbols(sampled_values >= -4/3 & sampled_values < 0) = 1;
received_symbols(sampled_values <= -4/3)  = 0;

% Convert integers back to voltages safely for the time-domain plot at the end
received_voltages = -2 + (4/3) * received_symbols;
output_received_data = repelem(received_voltages, rect_length);

% make the the full continous waveform
%output_received_data = repelem(received_symbols, rect_length);

%% Powers & BER Rate 
signal_power = mean(signal_values.^2);
interf_power_sp = mean((PAM_sing_pole(:) - PAM(:)).^2);

SIR_sp_dB = 10 * log10(signal_power / interf_power_sp);
% after equalization
interference_power = residual_ISI_power;
SIR_dB = 10*log10(signal_power / interference_power);

num_errors = sum(received_symbols ~= symbol_values);
SER_simulated = num_errors / N_symbols;
BER_simulated(idx) = SER_simulated/2;
% Noise after complete receiver
enhanced_noise_power = noise_power*noise_enhancement_factor;
total_impairment = residual_ISI_power + enhanced_noise_power;
% residual isi should equal best mse 
SINR_linear_theoretical = signal_power * effective_gain^2/total_impairment;
%BER_theoretical(idx) = qfunc(sqrt(SINR_linear_theoretical));
% PAM-4 Theoretical BER Approximation: SER = 1.5 * Q(sqrt(0.2 * SINR)), BER = SER / 2
SER_theoretical = 1.5 * qfunc(sqrt(0.2 * SINR_linear_theoretical));
BER_theoretical(idx) = SER_theoretical/2;
% Print progress
fprintf('Eb/No = %2d dB | Errors = %4d | BER = %e | SINR = %2d dB \n', current_EbN0_dB, num_errors, BER_simulated(idx), SINR_linear_theoretical);

end
fprintf('--- Pre-Channel ---------\n');
fprintf('Signal Power = %.4f\n\n', signal_power);
fprintf('--- Pre-Equalization SIR ---------\n\n');
fprintf('Interference Power = %.6f\n', interference_power);
fprintf('Single Pole Channel SIR = %.2f dB\n\n', SIR_sp_dB);
fprintf('--- AFTER-Equalization SIR ---------\n\n');
fprintf('Signal-to-Interference Ratio (SIR) = %.2f dB\n', SIR_dB);

fprintf('--- BER Results ---------\n');
fprintf('Total Bit Errors   = %d out of %d\n', num_errors, N_symbols);
fprintf('Simulated BER      = %e\n', BER_simulated);
fprintf('Theoretical BER    = %e\n', BER_theoretical);

%% Plots 
figure;
stairs(time1, PAM, 'Color', "g", 'LineStyle','--', 'LineWidth', 6.5 ,'DisplayName','no channel');
title('timedomain Comparison' );
hold on;
plot(time1, PAM_sing_pole, 'Color', "b", 'LineWidth', 2.5, 'DisplayName','Single pole channel');
hold on;
plot(time1, saved_very_noisy_equalized, 'Color', "r", 'LineWidth', 1.5, 'DisplayName','Very Noisy Equalized (Eb/No = 0 dB)');
hold on;
plot(time1, equalized_output_filtered, 'Color', "c", 'LineWidth', 2.5, 'DisplayName','equalized and filtered effect');

xlim([0 Tb*10]);
ylim([-5 5]);
grid on;
legend;

figure;
semilogy(EbN0_dB, BER_theoretical, 'r--', 'LineWidth', 2);
hold on;
semilogy(EbN0_dB, BER_simulated, 'bo', 'LineWidth', 1.5, 'MarkerSize', 8);
grid on;
title('BER vs E_b/N_0');
xlabel('E_b/N_0 (dB)');
ylabel('Bit Error Rate (BER)');
legend('Theoretical (Q-function)', 'Simulated (Monte Carlo)');
axis([min(EbN0_dB) max(EbN0_dB) 1e-5 1]);


figure;
plot(freq_domain, abs(single_pole_out));
hold on 
plot(freq_domain, low_pass_filter_fd,'LineWidth', 10 );
title("Low pass filter effect on traiing");
grid on;


%% Plotting Magnitude Responses
figure;
plot(freq_domain, 20*log10(abs(tran_fun)), 'r--', 'LineWidth', 2, 'DisplayName', 'Low-Pass Channel (tran\_fun)');
hold on;
plot(freq_domain, 20*log10(abs(equalizer_fd)), 'b-', 'LineWidth', 2, 'DisplayName', 'High-Pass Equalizer (equalizer\_fd)');

grid on;
title('Frequency Response: Single-Pole Channel vs. LS Equalizer');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
legend('Location', 'southwest');

% Limit x-axis to positive frequencies up to Nyquist (Fs/2)
xlim([0 main_lobe_freq/2]);
ylim([-30 15]);

%% eye diagram
eyediagram(PAM , 2 * rect_length);
eyediagram(saved_very_noisy_equalized2 , 2 * rect_length);
title('Eye Diagram: filtered Equalized Output ' );
% Plot Eye Diagram for the noisy 0 dB case
eyediagram(saved_very_noisy_equalized, 2*rect_length, 2*Tb, sampling_offset);
title('Eye Diagram of Very Noisy Equalized Signal (0 dB)');



