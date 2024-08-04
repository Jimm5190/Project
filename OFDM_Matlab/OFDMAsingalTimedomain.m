% Parameter
numSubcarriers = 8; % subcarrier
bandwidth = 20e6;   % bandwidth
subcarrierSpacing = bandwidth / numSubcarriers; % subcarrier space
fs = bandwidth;     % fs
nfft = 1024;        % FFT number
% subcarrier frequency
subcarrierFrequencies = (-numSubcarriers/2:numSubcarriers/2-1) * subcarrierSpacing;
% generate subcarrier
f = linspace(-fs/2, fs/2, nfft);
subcarriers = zeros(numSubcarriers, nfft);
for k = 1:numSubcarriers
    subcarriers(k, :) = abs(sinc((f - subcarrierFrequencies(k)) / subcarrierSpacing));
end
% plot
% Define a set of colors. You can modify or extend this list as needed.
colors = [
    0 0.4470 0.7410; % Blue
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
    0 0.4470 0.7410; 
];

% check
if numSubcarriers > size(colors, 1)
    error('Not enough colors defined for the number of subcarriers');
end

figure;
hold on;
for k = 1:numSubcarriers
    plot(f, subcarriers(k, :), 'Color', colors(k, :), 'DisplayName', ['Subcarrier ' num2str(k)]);
end
xlabel('Frequency (Hz)');
ylabel('Amplitude');
title('OFDM Subcarriers');
% legend show;
grid on;
hold off;

