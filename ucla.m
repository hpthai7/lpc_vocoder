%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EE214A - Digital Speech Processing - Class Project, Winter 2002
%
% Speech Coding using Linear Predictive Coding (LPC)
%
% Author: Philipp Steurer, 03/04/2002
%
% Project team: Daniel Thell, Ozgu Ozun, Philipp Steurer
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; % clear the command line
%clear all; % clear the workspace
close all;

%
% system constants
% ---------------
InputFilename = 'recorded_8000Hz';

[inspeech, Fs, bits] = wavread(InputFilename); % read the wavefile

outspeech1 = ucla_speechcoder(inspeech);
%outspeech2 = speechcoder2(inspeech);

% display the results
figure(1);
subplot(3,1,1);
plot(inspeech);
grid;
subplot(3,1,2);
plot(outspeech1);
grid;
% subplot(3,1,3);
% plot(outspeech2);
% grid;
% 
% disp('Press a key to play the original sound!');
% pause;
% soundsc(inspeech, 8000);
% 
% disp('Press a key to play the LPC compressed sound!');
% pause;
% soundsc(outspeech1, 8000);