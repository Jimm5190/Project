clear all;
close all;

% 設置參數
SNRindB1 = 0:2:15;    % 用於QAM的模擬SNR範圍
SNRindB2 = 0:0.1:15;  % 用於QAM的理論SNR範圍

M = 16;
k = log2(M);

% QAM模擬錯誤率
smld_err_prb_order_ber_qam = zeros(1, length(SNRindB1));
smld_err_prb_order_ser_qam = zeros(1, length(SNRindB1));
smld_err_prb_gray_ber_qam = zeros(1, length(SNRindB1));
smld_err_prb_gray_ser_qam = zeros(1, length(SNRindB1));
smld_err_prb_order_ber_psk = zeros(1, length(SNRindB1));
smld_err_prb_order_ser_psk = zeros(1, length(SNRindB1));
smld_err_prb_gray_ber_psk = zeros(1, length(SNRindB1));
smld_err_prb_gray_ser_psk = zeros(1, length(SNRindB1));
smld_err_prb_gray_ber_psk_rayleigh = zeros(1, length(SNRindB1));
smld_err_prb_gray_ser_psk_rayleigh = zeros(1, length(SNRindB1));
smld_err_prb_order_ber_apsk = zeros(1, length(SNRindB1));
%smld_err_prb_order_ser_apsk = zeros(1, length(SNRindB1));
smld_err_prb_gray_ber_apsk = zeros(1, length(SNRindB1));
%smld_err_prb_gray_ser_apsk = zeros(1, length(SNRindB1));
smld_err_prb_order_ber_apsk_comp = zeros(1, length(SNRindB1));
smld_err_prb_gray_ber_apsk_comp  = zeros(1, length(SNRindB1));

for i = 1:length(SNRindB1)

    N = 1000000;
    % 生成數據源
    dsource = randi([1 M], 1, N); % 生成 N 個 1 到 16 之間的隨機整數

    [ser, ber] = mapping(SNRindB1(i),'QAM','order',N,dsource); % 模擬錯誤率
    smld_err_prb_order_ser_qam(i) = ser; % 儲存符號錯誤率
    smld_err_prb_order_ber_qam(i) = ber; % 儲存比特錯誤率
    [ser, ber] = mapping(SNRindB1(i),'QAM','gray',N,dsource); 
    smld_err_prb_gray_ser_qam(i) = ser; 
    smld_err_prb_gray_ber_qam(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'PSK','order',N,dsource);
    smld_err_prb_order_ser_psk(i) = ser; 
    smld_err_prb_order_ber_psk(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'PSK','gray',N,dsource); 
    smld_err_prb_gray_ser_psk(i) = ser; 
    smld_err_prb_gray_ber_psk(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'APSK','order',N,dsource); 
    %smld_err_prb_order_ser_apsk(i) = ser; 
    smld_err_prb_order_ber_apsk(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'APSK','gray',N,dsource); 
    %smld_err_prb_gray_ser_apsk(i) = ser; 
    smld_err_prb_gray_ber_apsk(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'APSK1','order',N,dsource); 
    smld_err_prb_order_ber_apsk_comp(i) = ber; 
    [ser, ber] = mapping(SNRindB1(i),'APSK1','gray',N,dsource); 
    smld_err_prb_gray_ber_apsk_comp(i) = ber; 
    [ser, ber] = mapping_rayleigh(SNRindB1(i),'PSK','gray',N,dsource); 
    smld_err_prb_gray_ser_psk_rayleigh(i) = ser; 
    smld_err_prb_gray_ber_psk_rayleigh(i) = ber;
end

% 理論錯誤率
theo_err_prb_ser_qam = zeros(1, length(SNRindB2));
theo_err_prb_ber_qam = zeros(1, length(SNRindB2));
theo_err_prb_ser_gray_qam = zeros(1, length(SNRindB2));
theo_err_prb_ber_gray_qam = zeros(1, length(SNRindB2));
theo_err_prb_ser_psk = zeros(1, length(SNRindB2));
theo_err_prb_ber_psk = zeros(1, length(SNRindB2));
for i = 1:length(SNRindB2)
    SNR = 10^(SNRindB2(i) / 10); % 信噪比
    % QAM理論符號錯誤率
    theo_err_prb_ser_qam(i) = 4 * qfunc(sqrt(3 * k * SNR / (M - 1)));
    theo_err_prb_ber_qam(i) = qfunc(sqrt(3 * k * SNR / (M - 1)));
    theo_err_prb_ber_gray_qam(i) = 1/8*(2*erfc(sqrt(2*SNR/5))+erfc(3*sqrt(2*SNR/5))-erfc(5*sqrt(2*SNR/5)))+1/8*(erfc(sqrt(2*SNR/5))+erfc(3*sqrt(2*SNR/5)));
    theo_err_prb_ser_psk(i) = erfc(sqrt(k*10.^(SNRindB2(i)/10))*sin(pi/M));
    theo_err_prb_ber_psk(i) = (1/k)*erfc(sqrt(k*10.^(SNRindB2(i)/10))*sin(pi/M));
end

%figure qam gray and order
figure;
semilogy(SNRindB1, smld_err_prb_order_ser_qam, 's-');
hold on;
semilogy(SNRindB1, smld_err_prb_order_ber_qam, 's-');
semilogy(SNRindB1, smld_err_prb_gray_ser_qam, 'o-');
semilogy(SNRindB1, smld_err_prb_gray_ber_qam, 'o-');
semilogy(SNRindB2, theo_err_prb_ser_qam, '--');
semilogy(SNRindB2, theo_err_prb_ber_qam, '--');
semilogy(SNRindB2, theo_err_prb_ber_gray_qam, '--');
hold off;
legend('Simulated Order SER', 'Simulated Order BER', 'Simulated Gray SER','Simulated Gray BER','Theoretical SER', 'Theoretical BER','Theoretical Gray BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-QAM from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'qam.png';
print(fname, '-dpng');

%figure psk gray and order
figure
semilogy(SNRindB1, smld_err_prb_order_ser_psk, 's-');
hold on;
semilogy(SNRindB1, smld_err_prb_order_ber_psk, 's-');
semilogy(SNRindB1, smld_err_prb_gray_ser_psk, 'o-');
semilogy(SNRindB1, smld_err_prb_gray_ber_psk, 'o-');
semilogy(SNRindB2, theo_err_prb_ser_psk, '--');
semilogy(SNRindB2, theo_err_prb_ber_psk, '--');
hold off;
legend('Simulated Order SER', 'Simulated Order BER', 'Simulated Gray SER','Simulated Gray BER','Theoretical SER', 'Theoretical BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-PSK from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'PSK.png';
print(fname, '-dpng');

%compare qam and psk
figure;
semilogy(SNRindB1, smld_err_prb_order_ber_qam, 's-');
hold on;
semilogy(SNRindB1, smld_err_prb_gray_ber_qam, '^-');
semilogy(SNRindB1, smld_err_prb_order_ber_psk, 's-');
semilogy(SNRindB1, smld_err_prb_gray_ber_psk, '^-');
semilogy(SNRindB2, theo_err_prb_ber_qam, '--');
semilogy(SNRindB2, theo_err_prb_ber_gray_qam, '--');
semilogy(SNRindB2, theo_err_prb_ber_psk, '--');
hold off;
legend('Simulated Order QAM BER', 'Simulated Gray QAM BER','Simulated Order PSK BER','Simulated Gray PSK BER','Theoretical QAM BER','Theoretical QAM Gray BER', 'Theoretical PSK BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-QAM and 16-PSK BER from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'QAMPSK.png';
print(fname, '-dpng');

%compare qam、psk and apsk
figure;
semilogy(SNRindB1, smld_err_prb_order_ber_qam, 's-');
hold on;
semilogy(SNRindB1, smld_err_prb_gray_ber_qam, 's-');
semilogy(SNRindB1, smld_err_prb_order_ber_psk, 'o-');
semilogy(SNRindB1, smld_err_prb_gray_ber_psk, 'o-');
semilogy(SNRindB1, smld_err_prb_order_ber_apsk, '^-');
semilogy(SNRindB1, smld_err_prb_gray_ber_apsk, '^-');
hold off;
legend('Simulated Order QAM BER', 'Simulated Gray QAM BER','Simulated Order PSK BER','Simulated Gray PSK BER','Simulated Order APSK BER','Simulated Gray APSK BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-QAM、16-PSK and 16-APSK BER from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'QAMPSKAPSK.png';
print(fname, '-dpng')

%figure apsk angle
figure;
semilogy(SNRindB1, smld_err_prb_order_ber_apsk, '^-');
hold on;
semilogy(SNRindB1, smld_err_prb_gray_ber_apsk, '^-');
semilogy(SNRindB1, smld_err_prb_order_ber_apsk_comp, 'o-');
semilogy(SNRindB1, smld_err_prb_gray_ber_apsk_comp, 'o-');
hold off;
legend('Order APSK BER', 'Gray APSK BER','Order APSK Not rotation BER','Gray APSK Not rotation BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-APSK and Not rotation from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'APSK.png';
print(fname, '-dpng')

%Add Rayleigh fading
figure
semilogy(SNRindB1, smld_err_prb_gray_ser_psk, 's-');
hold on;
semilogy(SNRindB1, smld_err_prb_gray_ber_psk, 's-');
semilogy(SNRindB1, smld_err_prb_gray_ser_psk_rayleigh, 'o-');
semilogy(SNRindB1, smld_err_prb_gray_ber_psk_rayleigh, 'o-');
semilogy(SNRindB2, theo_err_prb_ser_psk, '--');
semilogy(SNRindB2, theo_err_prb_ber_psk, '--');
hold off;
legend('Simulated AWGN SER', 'Simulated AWGN BER', 'Simulated Rayleigh SER','Simulated Rayleigh BER','Theoretical SER', 'Theoretical BER');
xlabel('E_b/N_0 in dB', 'fontsize', 16, 'fontname', 'Helvetica');
ylabel('Error Probability', 'fontsize', 16, 'fontname', 'Helvetica');
title('Performance of 16-PSK by Gray code from Monte Carlo simulation', 'fontsize', 12, 'fontname', 'Helvetica');
fname = 'Rayleigh.png';
print(fname, '-dpng');

function [ser, ber] = mapping(snr_in_dB,mod_type,mapping_name,N,dsource)
% [ser, ber] = qam(snr_in_dB)
%       計算給定 snr_in_dB 下的符號錯誤率 (SER) 和比特錯誤率 (BER)

d = 1;                      % 符號之間的最小距離
Eav = 10 * d^2;             % 每個符號的能量
snr = 10^(snr_in_dB / 10);  % SNR per bit (給定)
sgma = sqrt(Eav / (8 * snr)); % 噪聲的標準差
M = 16;
r = 3.5;
% 映射到信號星座圖
if strcmp(mod_type,'QAM')
    mapping = [-3*d 3*d; -d 3*d; d 3*d; 3*d 3*d;
               -3*d d; -d d; d d; 3*d d;
               -3*d -d; -d -d; d -d; 3*d -d;
               -3*d -3*d; -d -3*d; d -3*d; 3*d -3*d];
elseif strcmp(mod_type,'PSK') % 16PSK
    mapping = zeros(M, 2);
        for i = 0:M-1
            theta = 2 * pi * i / M;
            mapping(i+1, :) = [cos(theta), sin(theta)] * r*d;
        end
elseif strcmp(mod_type,'APSK') %APSK
    mapping = zeros(M, 2);
        for i = 0:3
            theta = (2 * i + 1) * pi / 4 ;
            mapping(i+1, :) = [cos(theta), sin(theta)] * 1*d;
        end
        for j = 4:M-1;
            theta = 2 * pi * (j-4) / (M-4) + pi / (M-4);
            mapping(j+1, :) = [cos(theta), sin(theta)] * r*d;
        end   
else %APSK offset
        mapping = zeros(M, 2);
        for i = 0:3
            theta = (2 * i + 1) * pi / 4 ;
            mapping(i+1, :) = [cos(theta), sin(theta)] * 1*d;
        end
        for j = 4:M-1;
            theta = 2 * pi * (j-4) / (M-4) ;
            mapping(j+1, :) = [cos(theta), sin(theta)] * r*d;
        end
end

% Gray 編碼對應的比特映射
if strcmp(mapping_name,'gray')
    if strcmp(mod_type,'QAM')
        code_mapping = [1 0 1 1; 1 0 0 1; 0 0 0 1; 0 0 1 1;
                     1 0 1 0; 1 0 0 0; 0 0 0 0; 0 0 1 0;
                     1 1 1 0; 1 1 0 0; 0 1 0 0; 0 1 1 0;
                     1 1 1 1; 1 1 0 1; 0 1 0 1; 0 1 1 1];
    elseif strcmp(mod_type,'PSK')
        code_mapping = [0 0 0 0; 0 0 0 1; 0 0 1 1; 0 0 1 0;
                     0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0;
                     1 1 0 0; 1 1 0 1; 1 1 1 1; 1 1 1 0;
                     1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];
    else %APSK
        code_mapping = [1 1 0 0; 1 1 1 0; 1 1 1 1; 1 1 0 1;
                     0 1 0 0; 0 0 0 0; 1 0 0 0; 1 0 1 0;
                     0 0 1 0; 0 1 1 0; 0 1 1 1; 0 0 1 1;
                     1 0 1 1; 1 0 0 1; 0 0 0 1; 0 1 0 1];
    end
    else
         code_mapping = [0 0 0 0; 0 0 0 1; 0 0 1 0; 0 0 1 1;
                        0 1 0 0; 0 1 0 1; 0 1 1 0; 0 1 1 1;
                        1 0 0 0; 1 0 0 1; 1 0 1 0; 1 0 1 1;
                        1 1 0 0; 1 1 0 1; 1 1 1 0; 1 1 1 1];
end
qam_sig = mapping(dsource, :);

% 接收信號
n = sgma * randn(N, 2); % 生成高斯噪聲
r = qam_sig + n;        % 加上噪聲生成接收信號

% 檢測和錯誤概率計算
num_of_symbol_errors = 0;
num_of_bit_errors = 0;

for i = 1:N
    min_distance = inf; % 初始化最小距離為無限大
    for k = 1:M
        % 計算接收信號和每個星座點之間的歐幾里得距離
        distance = norm(r(i, :) - mapping(k, :));
        if distance < min_distance
            decis = k; % 如果當前距離小於之前的最小距離，則更新決策
            min_distance = distance;
        end
    end
    if decis ~= dsource(i)
        num_of_symbol_errors = num_of_symbol_errors + 1;
        
        % 計算原始符號和決策符號之間的比特錯誤數
        original_bits = code_mapping(dsource(i), :);
        decided_bits = code_mapping(decis, :);
        bit_errors = sum(original_bits ~= decided_bits);
        num_of_bit_errors = num_of_bit_errors + bit_errors;
    end
end

ser = num_of_symbol_errors / N;  % 計算符號錯誤率
ber = num_of_bit_errors / (N * log2(M)); % 計算比特錯誤率

end

function [ser, ber] = mapping_rayleigh(snr_in_dB, mod_type, mapping_name, N, dsource)
    % 計算給定 snr_in_dB 下的符號錯誤率 (SER) 和比特錯誤率 (BER)
    d = 1;                      % 符號之間的最小距離
    Eav = 10 * d^2;             % 每個符號的能量
    snr = 10^(snr_in_dB / 10);  % SNR per bit (給定)
    sgma = sqrt(Eav / (8 * snr)); % 噪聲的標準差
    M = 16;
    r = 3.5;

    % 映射到信號星座圖
    mapping = zeros(M, 2);
    for i = 0:M-1
        theta = 2 * pi * i / M;
        mapping(i+1, :) = [cos(theta), sin(theta)] * r*d;
    end
    % Gray 編碼對應的比特映射
    code_mapping = [0 0 0 0; 0 0 0 1; 0 0 1 1; 0 0 1 0;
                    0 1 1 0; 0 1 1 1; 0 1 0 1; 0 1 0 0;
                    1 1 0 0; 1 1 0 1; 1 1 1 1; 1 1 1 0;
                    1 0 1 0; 1 0 1 1; 1 0 0 1; 1 0 0 0];

    qam_sig = mapping(dsource, :);

    % Rayleigh衰落
    sigma = 1;  % Rayleigh分佈的參數
    alpha = sigma * sqrt(-2 * log(rand(N, 1)));  % 生成Rayleigh隨機變數

    % 接收信號
    n = sgma * randn(N, 2); % 生成高斯噪聲
    r = (qam_sig .* alpha) + n; % 加上噪聲生成接收信號

    % 檢測和錯誤概率計算
    num_of_symbol_errors = 0;
    num_of_bit_errors = 0;

    for i = 1:N
        min_distance = inf; % 初始化最小距離為無限大
        for k = 1:M
            % 計算接收信號和每個星座點之間的歐幾里得距離
            distance = norm(r(i, :) - mapping(k, :));
            if distance < min_distance
                decis = k; % 如果當前距離小於之前的最小距離，則更新決策
                min_distance = distance;
            end
        end
        if decis ~= dsource(i)
            num_of_symbol_errors = num_of_symbol_errors + 1;
            
            % 計算原始符號和決策符號之間的比特錯誤數
            original_bits = code_mapping(dsource(i), :);
            decided_bits = code_mapping(decis, :);
            bit_errors = sum(original_bits ~= decided_bits);
            num_of_bit_errors = num_of_bit_errors + bit_errors;
        end
    end

    ser = num_of_symbol_errors / N;  % 計算符號錯誤率
    ber = num_of_bit_errors / (N * log2(M)); % 計算比特錯誤率
end

