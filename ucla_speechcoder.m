function [ outspeech ] = ucla_speechcoder( inspeech )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% EE214A - Digital Speech Processing - Class Project, Winter 2002
%
% Speech Coding using Linear Predictive Coding (LPC)
% The desired order can be selected in the system constants section.
% For the excitation impulse-trains are used. The result does not sound very
% well but with this solution it is possible to achieve a low bitrate!
%
% Author: Philipp Steurer, 03/04/2002
%
% Parameters:
% inspeech : wave data with sampling rate Fs
% (Fs can be changed underneath if necessary)
%
% Returns:
% outspeech : wave data with sampling rate Fs
% (coded and resynthesized)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% arguments check
% ---------------
if ( nargin ~= 1)
    error('argument check failed');
end;

%
% system constants
% ----------------
Fs = 8000; % sampling rate in Hertz (Hz)
Order = 18; % order of the model used by LPC

%
% main
% ----

% encoded the speech using LPC
[aCoeff, resid, pitch, G, parcor, stream] = ucla_proclpc(inspeech, Fs, Order);

% decode/synthesize speech using LPC and impulse-trains as excitation
outspeech = ucla_synlpc(aCoeff, pitch, Fs, G);