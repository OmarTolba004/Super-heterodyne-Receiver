%{
 *  FILE DESCRIPTION
 *  -------------------------------------------------------------------------------------------------------------------
 *  File:  		  SuperHeterodyne.m
 *
 *  Description:  MATLAB file for Super-heterodyne Receiver Project
 *
 *  -------------------------------------------------------------------------------------------------------------------
 *	Author: 	  Omar Tolba & Omar Mustafa
 *	Date:		  15/12/2022
%}
%% Reading Audio Signal
% clear all WorkSpace Variables and Command Window
clear; clc;
% Reading the audio file and storing its data
[audioSample1,samplingFrequency1] = audioread('Short_WRNArabic.wav');
[audioSample2,samplingFrequency2] = audioread('Short_SkyNewsArabia.wav');
% converting the two channels into monophopic
audioSample1 = audioSample1(:,1)+audioSample1(:,2);
audioSample2 = audioSample2(:,1)+audioSample2(:,2);


%% Calculating Base band BW
% Specifying the length of FFT
N=2^20;
% Applying FFT
Y1=fft(audioSample1,N);
Y2=fft(audioSample2,N);
% Get the positive and negative frequencies
k=-N/2:N/2-1;
% Map it to actual frequencies --> note (samplingFrequency1==samplingFrequency2)
z=k*samplingFrequency1/N;
% Plotting FFT output against actual frequecnies
figure(1);
subplot(2,1,1);
plot(z,fftshift(abs(Y1)))
title('First Audio Signal After Applying FFT'); xlabel('Frequecny in Hz');
subplot(2,1,2);
plot(z,fftshift(abs(Y2)))
title('Second Audio Signal After Applying FFT'); xlabel('Frequecny in Hz');
% We can See that Audio Signal BaseBand Bw = 22 KHz

%% Resampling  
% Sampling interval --> note (samplingFrequency1==samplingFrequency2)
Ts = 1/samplingFrequency1;

% Resample audio data at a higher rate
resamplingFactor = 10;
audioSample1_ = interp(audioSample1,resamplingFactor);
audioSample2_ = interp(audioSample2,resamplingFactor);

% Time Intervals --> note : length(audioSample1)*resamplingFactor == length(audioSample1_)
% and length(audioSample2)*resamplingFactor == length(audioSample2_)
% and t1!= t2
t1 = linspace(0,Ts*length(audioSample1_),length(audioSample1)*resamplingFactor);
t2 = linspace(0,Ts*length(audioSample2_),length(audioSample2)*resamplingFactor);

% Plotting audio samples after resampling
figure(2);
subplot(2,1,1)
plot(t1,audioSample1_);
title('Audio Samples for first signal after resampling');
subplot(2,1,2)
plot(t2,audioSample2_);
title('Audio Samples for second signal after resampling');

%% AM Modulatoion
% Carrier Frequecny for first signal
fc1 = 10e4;
% Carrier Frequecny for second signal
fc2 = 10e4+5e4;
% Carrier for first Signal
yc1 = cos(2*pi*fc1*t1);
% Carrier for second Signal
yc2 = cos(2*pi*fc2*t2);
% Modulation for first Signal
sm1 = yc1.*audioSample1_';
% Modulation for Second Signal
sm2 = yc2.*audioSample2_';
% Analyzing in Freq. Spectrum for first audio singal
sa0=dsp.SpectrumAnalyzer('SampleRate',samplingFrequency1*10);
step(sa0,sm1');
% Analyzing in Freq. Spectrum for second audio singal
sa1=dsp.SpectrumAnalyzer('SampleRate',samplingFrequency2*10);
step(sa1,sm2');

%% padding the short signals with zeros so they have all equal length.
m1=length(sm1);
m2=length(sm2);
if m1>=m2
    for i=m2:m1
        sm2(1,i)=0;
    end
else 
     for i=m1:m2
         sm1(1,i)=0;
     end
end

%% Multiplexing the two audio signals
sm_t = sm1+sm2;
% Analyzing in Freq. Spectrum for multiplixed signal
sa3=dsp.SpectrumAnalyzer('SampleRate',samplingFrequency2*10);
step(sa3,sm_t');

