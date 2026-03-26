clear all;
close all;

%Step 1: Define Simulation Parameters
kvelocityoflight = 3*1e8; % (m/s)
carrier_freq_hz = 2*1e9;
waveleigh = kvelocityoflight / carrier_freq_hz;
velocity = 30; % m/s
max_fd = velocity / waveleigh;
sample_freq_hz = 20 * max_fd; % fs >= 10 max_fd
time = 5; % s
total_number_of_samples = time * sample_freq_hz; % Ns

%Step 2: Generate Gaussian Noise with Hermitian Symmetry
[W1, W2] = generate_complex_gaussian_noise(total_number_of_samples);

%Step 3: Apply the Doppler Filter
[W1_filter, W2_filter] = doppler_filter(total_number_of_samples, max_fd, W1, W2, sample_freq_hz);

%Step 4: Transform to the Time Domain
h1 = ifft(W1_filter);
h2 = ifft(W2_filter);

%Step 5: Form the Rayleigh Fading Envelope
r_n = h1+ 1i * h2;
r_n_envelope = sqrt(h1.^2 + h2.^2);
r_n_dB = 20 * log10(r_n_envelope);


%theta_n = atan2(h2, h1) * (180 / pi); 

%Step 6: Analyze the Simulated Channel
%Get the rayleigh pdf and compare with the theoretical
mean_real = mean(real(r_n));
var_real = var(real(r_n));

mean_imag = mean(imag(r_n));
var_imag = var(imag(r_n));
disp(['Mean of real part: ', num2str(mean_real)]);
disp(['Variance of real part: ', num2str(var_real)]);
disp(['Mean of imaginary part: ', num2str(mean_imag)]);
disp(['Variance of imaginary part: ', num2str(var_imag)])

num_bins = 50; 
[counts, bin_centers] = hist(r_n_envelope, num_bins); 
bin_width = bin_centers(2) - bin_centers(1); 
simulated_pdf = counts / (sum(counts) * bin_width); 

sigma = sqrt(mean(r_n_envelope.^2) / 2); % 計算 sigma
theoretical_pdf = (bin_centers / sigma^2) .* exp(-bin_centers.^2 / (2 * sigma^2));

rayleigh_mean = sigma * sqrt(pi / 2);
rayleigh_var = ((4 - pi) * sigma^2) / 2;
mean_sim = mean(r_n_envelope);
var_sim = var(r_n_envelope);
disp(['Simulation Rayleigh Mean: ', num2str(mean_sim)]);
disp(['Simulation Rayleigh Variance: ', num2str(var_sim)]);
disp(['Theoretical Rayleigh Mean: ', num2str(rayleigh_mean)]);
disp(['Theoretical Rayleigh Variance: ', num2str(rayleigh_var)]);

%Estimate the PSD of r(n)
window = hamming(1024); 
noverlap = 512; 
nfft = 2048; 
[pxx, f] = pwelch(r_n, window, noverlap, nfft, sample_freq_hz);
pxx = fftshift(pxx);  
f = (-nfft/2:nfft/2-1) * (sample_freq_hz / nfft); 

%autocorrelation
max_lag = floor(length(real(r_n)) / 2); 
[acf, lags] = autocorr(real(r_n),'NumLags',max_lag); 
dt = 1 / sample_freq_hz; 
time_lags = lags * dt; 
theoretical_acf = besselj(0, 2 * pi * max_fd * time_lags); 


%Generate Frequency-Domain Noise with Hermitian Symmetry function
function [W1, W2] = generate_complex_gaussian_noise(Ns)
    W1_raw = (randn(Ns, 1) + 1i * randn(Ns, 1)) / sqrt(2);
    W2_raw = (randn(Ns, 1) + 1i * randn(Ns, 1)) / sqrt(2);

    W1 = zeros(Ns, 1);
    W2 = zeros(Ns, 1);

    for k = 1:Ns
        if (k == 1) || (mod(Ns, 2) == 0 && k == Ns/2 + 1) 
            W1(k) = real(W1_raw(k));
            W2(k) = real(W2_raw(k));
        elseif k <= Ns/2
            W1(k) = W1_raw(k);
            W2(k) = W2_raw(k);
        else
            W1(k) = conj(W1_raw(Ns - k + 2)); 
            W2(k) = conj(W2_raw(Ns - k + 2)); 
        end
    end
end

% Doppler Spectrum Shaping function
function [W1_filter, W2_filter] = doppler_filter(Ns, max_fd, W1, W2, f_s)
    fk = (-Ns/2:Ns/2-1) * f_s / Ns;
    S_f = zeros(Ns, 1);
    
    for k = 1:Ns
        if abs(fk(k)) <= max_fd
            S_f(k) = 1 / (pi * max_fd * sqrt(max(1 - (fk(k)/max_fd)^2, 1e-3))); 
        else
            S_f(k) = 0;
        end
    end
    
    H_f = sqrt(S_f); 
    
    W1_filter = W1 .* H_f;
    W2_filter = W2 .* H_f;
end


% plot autocorrelation and bessel function
figure;
plot(time_lags, acf, 'b-', 'DisplayName', 'Computed Autocorrelation');
hold on;
plot(time_lags, theoretical_acf, 'r--', 'DisplayName', 'Jakes'' Autocorrelation');
xlabel('Time Delay (s)');
ylabel('Autocorrelation');
title('Comparison of Computed Autocorrelation and Theoretical Jakes'' Autocorrelation');
legend;
grid on;
xlim([0, 0.1]); 

% plot PSD
figure;
plot(f, 10*log10(pxx)); % 用dB顯示PSD
xlabel('Frequency (Hz)');
ylabel('Power Spectral Density (dB/Hz)');
title('Estimated Power Spectral Density of Complex r(n)');
xlim([-2500 2500]); 
grid on;


% plot h1
figure;
subplot(3,1,1);
plot(h1, 'b');
title('h1 - Complex Rayleigh Fading Signal');
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(real(h1), 'r');
title('Real Part of h1');
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(imag(h1), 'g');
title('Imaginary Part of h1');
xlabel('Samples');
ylabel('Amplitude');
grid on;

% plot h2
figure;
subplot(3,1,1);
plot(h2, 'b');
title('h2 - Complex Rayleigh Fading Signal');
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
plot(real(h2), 'r');
title('Real Part of h2');
xlabel('Samples');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(imag(h2), 'g');
title('Imaginary Part of h2');
xlabel('Samples');
ylabel('Amplitude');
grid on;

%plot |r(n)|
figure;
plot(1:total_number_of_samples, r_n_dB, 'b');
title('Correlated Rayleigh Magnitude Response', 'FontSize', 14);
xlabel('Samples');
ylabel('Amplitude (dB)');
grid on;

% plot rayleigh pdf
figure;
plot(bin_centers, simulated_pdf, 'bo-', 'DisplayName', 'Simulated PDF'); 
hold on;
plot(bin_centers, theoretical_pdf, 'r--', 'LineWidth', 2, 'DisplayName', 'Theoretical Rayleigh PDF'); 
xlabel('Magnitude');
ylabel('Probability Density');
title('Comparison of Simulated and Theoretical Rayleigh PDF');
legend;
grid on;

%{
figure;
plot(1:total_number_of_samples, theta_n, 'b');
title('Phase Response of Rayleigh Fading Signal', 'FontSize', 14);
xlabel('Samples');
ylabel('Phase (Degrees)');
ylim([-180, 180]); 
grid on;
%}