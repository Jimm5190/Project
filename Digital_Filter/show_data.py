import os
import numpy as np
import matplotlib.pyplot as plt


# 定義函數，繪製脈衝響應圖
def plot_impulse_response(file_path):
    # 讀取脈衝響應數據
    with open(file_path, 'r') as file:
        impulse_response = np.array([float(value) for value in file.read().splitlines()])

    # 生成時間軸
    time = np.arange(0, len(impulse_response))
    color = 'b'

    # 繪製脈衝響應圖
    plt.stem(time, impulse_response, linefmt=f'{color}-', basefmt=f'{color}-')

    plt.xlabel('Time')
    plt.ylabel('Amplitude')
    plt.title(f'Impulse Response - {os.path.basename(file_path).replace(".txt", "")}')

    # 將圖片保存為 JPEG 文件
    plt.savefig(f"{os.path.basename(file_path).replace('.txt', '')}.jpg")
    plt.close()


# 定義函數，繪製對數刻度頻譜圖
def plot_log_spectrum(file_path):
    sampling_rate = 48000  # 取樣频率為48000
    num_samples = 1202
    # 讀取頻譜數據
    with open(file_path, 'r') as file:
        # 將數據轉換為NumPy數組
        data = np.array([float(value) for value in file.read().splitlines()])[:num_samples//2]

    # 生成频率轴
    frequency_axis = np.fft.fftfreq(num_samples, d=1/sampling_rate)[:num_samples//2]

    # 繪製 log spectrum
    plt.plot(frequency_axis, data)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Magnitude log (dB)')
    
    # 標記頻率 1500 和 3500 的位置
    plt.axvline(x=1500, color='0.8', linestyle='--', label='1500 Hz')
    plt.axvline(x=3500, color='0.8', linestyle='--', label='3500 Hz')
    
    plt.title(f'Log Spectrum - {os.path.basename(file_path).replace(".txt", "")}')
    plt.legend()  # 顯示標記的圖例
    plt.savefig(f"{os.path.basename(file_path).replace('.txt', '')}.jpg")
    plt.close()



file_hL=("hL-8.txt", "hL-32.txt", "hL-1024.txt")
file_hR=("hR-8.txt", "hR-32.txt", "hR-1024.txt")
file_YL=("YL-8.txt", "YL-32.txt", "YL-1024.txt")
file_YR=("YR-8.txt", "YR-32.txt", "YR-1024.txt")

for i in range (3):
    plot_impulse_response(file_hL[i])
    plot_impulse_response(file_hR[i])
    plot_log_spectrum(file_YL[i])
    plot_log_spectrum(file_YR[i])
