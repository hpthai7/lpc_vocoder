clc
close all
clear all

% Record voice at sampling frequency of 8KHz in 5s
Fs = 4000;
duration = 5;
[y, fname] = record_voice(Fs, duration);

% Read
[y, Fsr] = wavread(fname);
wavplay(y, Fsr);
filename_mat = sprintf('recorded_%dHz.mat\n', Fsr);
save y filename_mat;

% Display
plot(y);
grid on;
legend('Signal sampled at 8KHz');
xlabel('Time (sec)');
ylabel('Amplitude');