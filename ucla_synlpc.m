function synWave = ucla_synlpc(aCoeff,pitch,sr,G,fr,fs,preemp)
% USAGE: synWave = synlpc(aCoeff,pitch,sr,G,fr,fs,preemp);
%
% This function synthesizes a (speech) signal based on a LPC (linear-
% predictive coding) model of the signal. The LPC coefficients are a
% short-time measure of the speech signal which describe the signal as the
% output of an all-pole filter. This all-pole filter provides a good
% description of the speech articulators; thus LPC analysis is often used
% in speech recognition and speech coding systems. The LPC analysis is done
% using the proclpc routine. This routine can be used to verify that the
% LPC analysis produces the correct answer, or as a synthesis stage after
% first modifying the LPC model.
%
% The results of LPC analysis are a new representation of the signal s(n) =
% G e(n) - sum from 1 to L a(i)s(n-i) where s(n) is the original data. a(i)
% and e(n) are the outputs of the LPC analysis with a(i) representing the
% LPC model. The e(n) term represents either the speech source's
% excitation, or the residual: the details of the signal that are not
% captured by the LPC coefficients. The G factor is a gain term.
%
% LPC synthesis produces a monaural sound vector (synWave) which is sampled
% at a sampling rate of "sr". The following parameters are mandatory aCoeff
% - The LPC analysis results, a(i). One column of L+1 numbers for each
% frame of data. The number of rows of aCoeff determines L. G - The LPC
% gain for each frame. pitch - A frame-by-frame estimate of the pitch of
% the signal, calculated by finding the peak in the residual's
% autocorrelation for each frame.
%
% The following parameters are optional and default to the indicated
% values. fr - Frame time increment, in ms. The LPC analysis is done
% starting every fr ms in time. Defaults to 20ms (50 LPC vectors a second)
% fs - Frame size in ms. The LPC analysis is done by windowing the speech
% data with a rectangular window that is fs ms long. Defaults to 30ms
% preemp - This variable is the epsilon in a digital one-zero filter which
% serves to preemphasize the speech signal and compensate for the 6dB per
% octave rolloff in the radiation function. Defaults to .9378.
%
% This code was graciously provided by: Delores Etter (University of
% Colorado, Boulder) and Professor Geoffrey Orsak (Southern Methodist
% University) It was first published in Orsak, G.C. et al. "Collaborative
% SP education using the Internet and MATLAB" IEEE SIGNAL PROCESSING
% MAGAZINE Nov. 1995. vol.12, no.6, pp. 23-32. Modifications by Philipp
% Steurer: Using impulse-trains for the voice excitation.

% (c) 1998 Interval Research Corporation A more complete set of routines
% for LPC analysis can be found at
% http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html

if (nargin < 5), fr = 20; end;
if (nargin < 6), fs = 30; end;
if (nargin < 7), preemp = .9378; end;

msfs = round(sr*fs/1000); % framesize in samples
msfr = round(sr*fr/1000); % framerate in samples
msoverlap = msfs - msfr;
ramp = [0:1/(msoverlap-1):1]';
[L1 nframe] = size(aCoeff); % L1 = 1+number of LPC coeffs

for frameIndex=1:nframe
    A = aCoeff(:,frameIndex);
    % first check if it is voiced or unvoiced sound:
    if ( pitch(frameIndex) ~= 0 )
        t = 0 : 1/sr : fs*10^(-3); % sr sample freq. for fr ms
        d = 0 : 1/pitch(frameIndex) : 1; % 1/pitchfreq. repetition freq.
        residFrame = (pulstran(t, d, 'tripuls', 0.001))'; % sawtooth width of 0.001s
        residFrame = residFrame + 0.01*randn(msfs+1,1);
    else
        residFrame = [];
        for m = 1:msfs
            residFrame = [residFrame; randn];
        end % for
    end;
    
    synFrame = filter(G(frameIndex), A', residFrame); % synthesize speech from LPC coeffs
    if(frameIndex==1) % add synthesize frames using a trapezoidal window
        synWave = synFrame(1:msfr);
    else
        synWave = [synWave; overlap+synFrame(1:msoverlap).*ramp; ...
            synFrame(msoverlap+1:msfr)];
    end
    if(frameIndex==nframe)
        synWave = [synWave; synFrame(msfr+1:msfs)];
    else
        overlap = synFrame(msfr+1:msfs).*flipud(ramp);
    end
end;
% synWave = filter(1, [1 -preemp], synWave);