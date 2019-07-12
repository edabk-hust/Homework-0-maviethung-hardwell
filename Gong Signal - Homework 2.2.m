%% Prolog: specgong.m from Software Receiver Design text
filename = 'gong.wav' ; % name of wave file goes here
[ x , sr ] = audioread( filename ) ; % read in wavefile
Ts = 1/ sr ; siz=length ( x ) ; % sample interval and # of samples
N = 2^16; x = x(1:N)'; % length for analysis
sound(x , 1 / Ts ) % play sound , if sound card installed
time = Ts * ( 0 : length(x)-1); % establish time base for plotting
subplot (2, 1, 1) , plot ( time , x );
title ('Unfiltered signal in time domain'); % and plot top figure
magx = abs ( fft(x) ) ; % take FFT magnitude
ssf = (0:N/2-1)/(Ts*N) ; % establish freq base for plotting
subplot(2, 1, 2), plot(ssf, magx ( 1:N/2));
title ('Unfiltered signal in frequncy domain'); % plot mag spectrum
%% (a) filter design
%% Design a filter using firpm that will remove the two highest partials from this sound without affecting the lowest partial.
%% Nyquist frequency, i.e. half the sampling rate
fnyquist = sr/2;
%% Define the passband frequency in Hz
fpassband = 530;
%% Define the stopband frequency in Hz
fstopband = 585;
ctfrequencies = [0 fpassband fstopband fnyquist];
pmfrequencies = ctfrequencies / fnyquist;
%% Define the number of coefficients for the filter
filterlength = 100;
lowpassamps = [ 1 1 0 0 ];
lowpassfilter = firpm( filterlength, pmfrequencies, lowpassamps );
%% (b) filter signal
%% Use the filter command to process the gong.wav file with your filter.
y = filter(lowpassfilter, 1, x);
%% (c) plot the results
%% Take the FFT of the resulting signal and verify that the partial at 520 remains while the others are removed.
subplot(4, 1, 1); plot (time, y);
title ('Filtered signal in time domain');
magx1 = abs( fft( y ) );
subplot(2,1,2); plot (ssf, magx1( 1:N/2));
title ('Filtered signal in frequency domain');
sound (y, 1/Ts);
%% downsampling by 2
yDownSampled = y(1:2:length(y));
sound(yDownSampled, 1/Ts);
% Compute filtered spectrum and frequency vector as in specgong.m
magY = abs(fft(y));
N = length(magY);
f =(0:N/2-1)/(Ts*N);
subplot(2, 1, 1);  plot(f, magY(1:N/2));
xlabel('Frequency (Hz)');
title('Filtered Magnitude Spectrum');
% Compute downsampled spectrum
magYd = abs(fft(yDownSampled));
N2 = length(magYd);
f2 = (0:N2/2-1)/(Ts*N2);
subplot(2, 1, 2); plot(f2, magYd(1:N2/2));
xlabel('Frequency (Hz)');
title('Filtered & Downsampled Magnitude Spectrum');



