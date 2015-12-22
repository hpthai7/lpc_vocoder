function [pitch_list, r] = calculate_pitch(input, Fs)
% Input signal
% duration in seconds
% correlation array r
clc;
close all;
frame_length = 20 * 10^(-3); % second

% fenetre de n_frame_samples echantillon
n_total_samples = length(input)
duration = n_total_samples / Fs;
n_frame_samples = n_total_samples * frame_length / duration

n_windows = n_total_samples / n_frame_samples;

i = 0; j = 0;
lpc_order = 18; % order of the model used by LPC. lpc_order + 1 coefficients
correlations = zeros(n_windows, n_frame_samples);
pitch_list = ones(1, n_windows);
while (i+n_frame_samples <= n_total_samples)
    excerpt(1:n_frame_samples) = input(i+1:i+n_frame_samples);
    r = autocorr(excerpt);
    [peaks, locs] = findpeaks(r);
    i = i+n_frame_samples; j = j+1;
    
    correlations(j,:) = r;
    
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
% %        disp('Pas de pitch');
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
end

correlations_excerpt = correlations(:,1:lpc_order + 1);

% continuer faire de computation de levinson
% corr_size = size(correlations_excerpt);
% coefs = zeros(corr_size);
% for i = 1:corr_size(1)
%     lpc_coefs = [levinson_durbin(correlations_excerpt(i,:))];
%     coefs(i,:) = [1, lpc_coefs];
% end

coefs = zeros(n_windows, lpc_order + 1);
for i = 1:n_windows
    % estimation error
    window = input((i-1)*n_frame_samples+1:i*n_frame_samples);
    errSig = filter(coefs(i,:), 1, window);

    % Levinson's method
    autoCorVec = correlations(i,1:lpc_order+1);
    err(1) = autoCorVec(1);
    k(1) = 0;
    A = [];
    for index=1:lpc_order
        numerator = [1 A.']*autoCorVec(index+1:-1:2)';
        denominator = -1*err(index);
        k(index) = numerator/denominator; % PARCOR coeffs
        A = [A+k(index)*flipud(A); k(index)];
        err(index+1) = (1-k(index)^2)*err(index);
    end
    
    coefs(i,:) = [1; A]';
end

% pitch detection % findpeak not effective
for i = 1:n_windows
    % estimation error
    window = input((i-1)*n_frame_samples+1:i*n_frame_samples);
    errSig = filter(coefs(i,:), 1, window);

%    G(nframe) = sqrt(err(L+1)); % gain
    autocorrErr = autocorr(errSig'); % calculate pitch & voicing information
    [B,I] = sort(autocorrErr);
    num = length(I);
    if B(num-1) > .01*B(num)
        pitch_list(i) = abs(I(num) - I(num-1));
        disp('a pitch');
    else
        pitch_list(i) = 0;
        disp('NOTTTTTT a pitch');
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
% SYTHESIS
% input: pitch_list (1 x windows), LPC coefficients (n_windows x n_window_size),

disp('Modelisation de source...');
% generate input for LPC filter
[n_windows n_coefs] = size(coefs);
n_excitations = n_frame_samples;
excitations = ones(n_windows, n_excitations)/10;
% for i = 1:n_windows
%     if pitch_list(i) == 0
%         white_noise = zeros(1, n_excitations);
%         excitations(i, :) = awgn(white_noise, 50); % snr
%     end
% end

disp('LPC (AR) filtering..');
Gain = 1;
% LPC (AR) filtering
idx = 0;
for i = 1:n_windows
%    coefs(i,1:n_coefs)
    ar_output = filter(Gain, -coefs(i,1:n_coefs), excitations(i,1:n_excitations));
    syn_output(idx + 1:idx + n_excitations, 1) = ar_output';
    idx = idx + n_excitations;
end

disp('Finish');
coef_size = size(coefs)
out_size = size(syn_output)

%wavplay(syn_output, Fs);

figure(1);
subplot(3,1,1);
plot(input);
grid;
subplot(3,1,2);
plot(syn_output);
grid;

disp('Press a key to play the original sound!');
pause;
soundsc(input, Fs);

disp('Press a key to play the LPC compressed sound!');
pause;
soundsc(syn_output, Fs);