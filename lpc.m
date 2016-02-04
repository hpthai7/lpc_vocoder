function lpc(input, Fs)
% Input signal
% Fs sampling frequency in Hz

clc;
close all;
frame_duration = 20 * 10^(-3); % second

% fenetre de n_frame_samples echantillon
n_total_samples = length(input)
duration = n_total_samples / Fs;
n_frame_samples = n_total_samples * frame_duration / duration

n_windows = n_total_samples / n_frame_samples;

i = 0; j = 0;
correlations = zeros(n_windows, n_frame_samples);
while (i+n_frame_samples <= n_total_samples)
    excerpt(1:n_frame_samples) = input(i+1:i+n_frame_samples);
    r = autocorr(excerpt);
    [peaks, locs] = findpeaks(r);
    i = i+n_frame_samples; j = j+1;
    
    correlations(j,:) = r;
end

% Estimate LPC coefficients using Levinson's method
lpc_order = 18; % order of LPC filter, we'll have lpc_order + 1 coefficients
coefs = zeros(n_windows, lpc_order + 1);
for i = 1:n_windows
    autocorrVec = correlations(i,1:lpc_order+1);
    err(1) = autocorrVec(1);
    k(1) = 0;
    A = [];
    for index=1:lpc_order
        numerator = [1 A.']*autocorrVec(index+1:-1:2)';
        denominator = -1*err(index);
        k(index) = numerator/denominator;
        A = [A+k(index)*flipud(A); k(index)];
        err(index+1) = (1-k(index)^2)*err(index);
    end
    coefs(i,:) = [1; A]';
end

% pitch detection % findpeak not effective
pitch_list = ones(1, n_windows);
for i = 1:n_windows
    % estimation error
    window = input((i-1)*n_frame_samples+1:i*n_frame_samples);
    error_signal = filter(coefs(i,:), 1, window);

    autocorrErr = autocorr(error_signal');
    [B,I] = sort(autocorrErr);
    num = length(I);
    if B(num-1) > .01*B(num)
        pitch_list(i) = abs(I(num) - I(num-1));
    else
        pitch_list(i) = 0;
    end
%     if pitch_list(i) ~= 0
%         [peaks, locs] = findpeaks(autocorrErr);
%         figure;
%         subplot(2,1,1); plot(errSig);legend('Windows');
%         subplot(2,1,2); plot(autocorrErr);legend('Correlation');
%         hold on;
%         plot(locs, peaks, 'r*');
%         hold off;
%     end
end

% find pitch % findpeak not effective
%     peaks = findpeaks(r);
%     if length(peaks) <= 2
%         pitch_list(j) = 0;
% %        disp('Pas de pitch');
%         continue;
%     end
%     sample_peaks = peaks(1:3);
%     sample_locs = [0, locs(1:3)];
%     pitches(1:3) = sample_locs(2:4) - sample_locs(1:3);
%     diff(1:2) = pitches(2:3) - pitches(1:2);
%     
%     max_difference = 8;
%     if abs(diff(1)) > max_difference || abs(diff(2)) > max_difference || diff(1)*diff(2) > 0
%         pitch_list(j) = 0;
%         disp('Pas de pitch');
%         continue;
%     end
    
%     figure;
%     subplot(2,1,1); plot(excerpt);legend('Windows');
%     subplot(2,1,2); plot(r);legend('Correlation');
%     hold on;
%     plot(locs, peaks, 'r*');
%     hold off;
%     sample_locs
%     pitches
%     diff


% SYTHESIS
% generate input for LPC filter
% input: pitch_list (1 x n_windows), coefs (n_windows x lpc_order),
excitations = zeros(n_windows, n_frame_samples+1);
for i = 1:n_windows
    if pitch_list(i) == 0
        white_noise = zeros(1, n_frame_samples+1);
        excitations(i, :) = awgn(white_noise, 50); % snr
    else
        t = 0 : 1/Fs : frame_duration; % 20ms
        d = 0 : 1/pitch_list(i) : 1; % 1/pitchfreq. repetition freq.
        pulses = (pulstran(t, d, 'tripuls', 0.001))'; % periodic excitation
        %white_noise = zeros(n_frame_samples+1, 1);
        %white_noise = awgn(white_noise, 50); % snr
        excitations(i, :) = pulses; % + white_noise;
        if (i == 1)
            plot(excitations(i,:));
        end
    end
end

% LPC (AR) filtering
Gain = 5;
idx = 0;
for i = 1:n_windows
%    coefs(i,1:n_coefs)
    ar_output = filter(Gain, coefs(i,:), excitations(i,:));
    size_ar = size(ar_output)
    syn_output(idx + 1:idx + n_frame_samples, 1) = ar_output(1:n_frame_samples)';
    idx = idx + n_frame_samples;
end

save 'synthesized_speech_8KHz.mat' syn_output;

figure(1);
subplot(2,1,1);
plot(input);
grid;
subplot(2,1,2);
plot(syn_output);
grid;

% FFT transform and frequency spectrum
NFFT = 2^nextpow2(n_total_samples); % Next power of 2 from length of y
Y = fft(input,NFFT)/n_total_samples;
f = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure;
plot(f,2*abs(Y(1:NFFT/2+1))) 
title('Single-sided amplitude spectrum of input signal')
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

Y2 = fft(syn_output,NFFT)/n_total_samples;

% Plot single-sided amplitude spectrum.
figure;
plot(f,2*abs(Y2(1:NFFT/2+1))) ;
title('Single-sided amplitude spectrum of synthesized signal');
xlabel('Frequency (Hz)');
ylabel('|Y(f)|');
grid on;

disp('Press a key to hear original speech');
pause;
soundsc(input, Fs);

disp('Press a key to hear synthesized speech');
pause;
soundsc(syn_output, Fs);