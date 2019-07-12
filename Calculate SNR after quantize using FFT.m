
close all; % closing all the open figures
clear all; % clearing all previous variables
%=
% Generate a simulation signal
%=
fin = 257;
fs = 8192;
N = 8192;
adc_resolution = 16; % Or Quantizer resolution

Nw = floor((fin*N)/fs);
amp = 0.5;
signal = amp*sin(2*pi*Nw/N*[0:N-1]);

figure(1);
subplot (2,1,1);
plot(signal);
title('Original Signal');
ylabel('Signal Amplitude');
xlabel('Time');

% Scalar quantization is implemented over here
no_quantiz_levels = power(2, adc_resolution);
quantize_signal = floor((no_quantiz_levels-1)*signal)/(no_quantiz_levels/2);%+no_quantiz_levels/2

subplot (2,1,2);
plot(signal-quantize_signal/2, 'g*');
title('Quantization error');

% FFT and SNR measurement
% Remember that hann window reduces the amplitude of the signal to half in
% frequency domain and that is why scaling is done by N/4 instead of N/2
window_type=ones(1,N); %hann(N)'; %Coherent sampling, no need for windowing
window_scaling_factor=sum(window_type);
fft_signal = fft(quantize_signal.*window_type,N)*(2/window_scaling_factor);
freq_scal = linspace (0, 0.5, N/2); % for normalized frequency
X=fft_signal.*conj(fft_signal);

figure(2);
plot(freq_scal, 10*log10(X(1:N/2)), '-x');
hold on
plot(freq_scal(Nw:Nw+2), 10*log10(abs(fft_signal(Nw:Nw+2))), 'r*');
title('FFT Plot (N = 8192, Hann window is used)');
xlabel('Frequency (Normalized)');

signal_value = (sum(fft_signal(Nw:Nw+2).*conj(fft_signal(Nw:Nw+2))));

noise_bins = [fft_signal(2:fin-1) fft_signal(fin+3:N/2)];
noise = (sum((abs(noise_bins)).^2)) / ((N/2)-3);
noise = 2*(sum(noise_bins.*conj(noise_bins)));
noise_formula = ((1/power(2,adc_resolution))^2)/12;

SNR = 10*log10(signal_value/noise)
SNR_byQuantizationFormula = 10*log10(signal_value/noise_formula)