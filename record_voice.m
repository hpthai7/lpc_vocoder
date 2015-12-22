function [y, str] = record_voice(Fs, duration)
% input: Hz, second
% output: signal y and filename str

disp(['RECORDING signal at ', num2str(Fs), 'Hz...']);
y = wavrecord(duration * Fs, Fs, 'int16');
wavplay(y, Fs);

str = sprintf('recorded_%dKHz.wav', Fs);
disp(['SAVING signal at ', num2str(Fs), 'KHz']);
wavwrite(y, Fs, 16, str); % 16bit
disp('FINISH');