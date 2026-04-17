clear all;
close all;

% Parameters
Nt = 2; % Transmit antennas
Nr_Alamouti = 1; % Receive antenna for Alamouti
Nr_Mux = 2; % Receive antennas for spatial multiplexing
M = 4; % QPSK
k = log2(M); % Bits per symbol
Nbits = 1e5;
EbNo_dB = 0:2:30;
EbNo = 10.^(EbNo_dB/10);
EsNo = EbNo * k;
No = 1 ./ EsNo;

% Results initialization
BER_Alamouti = zeros(1, length(EbNo_dB));
BER_ZF = zeros(1, length(EbNo_dB));
BER_MMSE = zeros(1, length(EbNo_dB));
Throughput_Alamouti = zeros(1, length(EbNo_dB));
Throughput_MMSE = zeros(1, length(EbNo_dB));
Throughput_ZF = zeros(1, length(EbNo_dB));

% QPSK symbol mapping table
symbol_map = [1+1j, -1+1j, 1-1j, -1-1j] / sqrt(2);
bit_map = [0 0; 0 1; 1 0; 1 1];

%% Alamouti STBC (2x1)
for i = 1:length(EbNo_dB)
    no_errors = 0;
    total_bits = 0;
    for iter = 1:Nbits/(2*k)
        % Generate 4 bits and map to 2 QPSK symbols
        bits = randi([0 1], 1, 4);
        idx1 = bi2de(bits(1:2), 'left-msb') + 1;
        idx2 = bi2de(bits(3:4), 'left-msb') + 1;
        s = [symbol_map(idx1), symbol_map(idx2)];

        % Alamouti encoding (2 time slots)
        s_Alamouti = [s(1), -conj(s(2)); s(2), conj(s(1))]; % [2 x 2]

        % Channel and noise (Nr = 1)
        h = (randn(1, Nt) + 1j*randn(1, Nt)) / sqrt(2); % [1 x 2]
        noise = sqrt(No(i)/2) * (randn(1,2) + 1j*randn(1,2)); % [1 x 2]

        % Transmit
        y = h * s_Alamouti + noise; % [1 x 2]

        % Alamouti combining
        y1 = y(1);
        y2 = y(2);
        r1 = conj(h(1)) * y1 + h(2) * conj(y2);
        r2 = conj(h(2)) * y1 - h(1) * conj(y2);
        r = [r1, r2];

        % Detection
        s_hat = zeros(1,2);
        for j = 1:2
            dists = abs(r(j) - symbol_map).^2;
            [~, idx] = min(dists);
            s_hat(j) = idx;
        end
        % Demap
        bits_hat = [bit_map(s_hat(1), :), bit_map(s_hat(2), :)];
        no_errors = no_errors + sum(bits_hat ~= bits);
        total_bits = total_bits + 4;
    end
    BER_Alamouti(i) = no_errors / total_bits;
    Throughput_Alamouti(i) = (k * 2 / 2) * (1 - BER_Alamouti(i)); % Rate-1/2
end

%% Spatial Multiplexing (2x2)
for i = 1:length(EbNo_dB)
    no_errors_ZF = 0;
    no_errors_MMSE = 0;
    total_bits = 0;
    for iter = 1:Nbits/(2*k)
        % Generate 4 bits and map to 2 QPSK symbols
        bits = randi([0 1], 1, 4);
        idx1 = bi2de(bits(1:2), 'left-msb') + 1;
        idx2 = bi2de(bits(3:4), 'left-msb') + 1;
        s = [symbol_map(idx1); symbol_map(idx2)];

        % Channel and noise
        H = (randn(Nr_Mux, Nt) + 1j*randn(Nr_Mux, Nt)) / sqrt(2);
        noise = sqrt(No(i)/2) * (randn(Nr_Mux,1) + 1j*randn(Nr_Mux,1));
        y = H * s + noise;

        % ZF detection
        s_ZF = pinv(H) * y;
        s_hat_ZF = zeros(1, 2);
        for j = 1:2
            dists = abs(s_ZF(j) - symbol_map).^2;
            [~, idx] = min(dists);
            s_hat_ZF(j) = idx;
        end
        bits_hat_ZF = [bit_map(s_hat_ZF(1), :), bit_map(s_hat_ZF(2), :)];
        no_errors_ZF = no_errors_ZF + sum(bits_hat_ZF ~= bits);

        % MMSE detection
        W = (H'*H + No(i)*eye(Nt)) \ H';
        s_MMSE = W * y;
        s_hat_MMSE = zeros(1, 2);
        for j = 1:2
            dists = abs(s_MMSE(j) - symbol_map).^2;
            [~, idx] = min(dists);
            s_hat_MMSE(j) = idx;
        end
        bits_hat_MMSE = [bit_map(s_hat_MMSE(1), :), bit_map(s_hat_MMSE(2), :)];
        no_errors_MMSE = no_errors_MMSE + sum(bits_hat_MMSE ~= bits);

        total_bits = total_bits + 4;
    end
    BER_ZF(i) = no_errors_ZF / total_bits;
    BER_MMSE(i) = no_errors_MMSE / total_bits;
    Throughput_MMSE(i) = 2 * k * (1 - BER_MMSE(i)); % Full rate
    Throughput_ZF(i) = 2 * k * (1 - BER_ZF(i));
end

%% Plot Results
figure;
semilogy(EbNo_dB, BER_Alamouti, '-o', 'LineWidth', 2); hold on;
semilogy(EbNo_dB, BER_ZF, '-*', 'LineWidth', 2);
semilogy(EbNo_dB, BER_MMSE, '-x', 'LineWidth', 2);
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('BER', 'FontSize', 12);
legend('Alamouti STBC (2x1)', 'ZF Detection (2x2)', 'MMSE Detection (2x2)', 'FontSize', 12);
title('BER Performance Comparison', 'FontSize', 14);
grid on;

figure;
plot(EbNo_dB, Throughput_Alamouti, '-o', 'LineWidth', 2); hold on;
plot(EbNo_dB, Throughput_ZF, '-x', 'LineWidth', 2);
plot(EbNo_dB, Throughput_MMSE, '-s', 'LineWidth', 2);
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Throughput (bits/symbol)', 'FontSize', 12);
legend('Alamouti', 'ZF ','MMSE', 'FontSize', 12);
title('Throughput vs. SNR', 'FontSize', 14);
grid on;

figure;
plot(EbNo_dB, Throughput_ZF, '-s', 'LineWidth', 2); hold on;
plot(EbNo_dB, Throughput_MMSE, '-s', 'LineWidth', 2);
xlabel('SNR (dB)', 'FontSize', 12);
ylabel('Throughput (bits/symbol)', 'FontSize', 12);
legend('ZF', 'MMSE ', 'FontSize', 12);
title('Throughput vs. SNR', 'FontSize', 14);
grid on;
