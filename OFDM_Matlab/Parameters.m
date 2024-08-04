%% System Parameters%%
SYS_Para.system_mode = 'OFDM'; % 'OFDM', 'OTFS', 'OTFS_LOW_COMP'
SYS_Para.carrier_freq_mhz = 2600;
SYS_Para.max_fft_size = 4096;
SYS_Para.max_subcarrier_spacing_khz = 480;
SYS_Para.unit_time_ms = 1 / (SYS_Para.max_fft_size * SYS_Para.max_subcarrier_spacing_khz); % tc
SYS_Para.cyclic_prefix_ratio = 9/128; % (CP length / OFDM Symbol Length)
SYS_Para.numerology.Bandwidth = readtable('5G_Numerology.xlsx','Sheet','Bandwidth');
SYS_Para.signal_bandwidth_mhz = 5; % 5,10,15,20,25,30,40,50,60,70,80,90,100
SYS_Para.subcarrier_spacing_khz = 30; % FR1 (BW<100MHz): 15,30,60 FR2 (BW<400MHz): 60,120,240
SYS_Para.subcarriers_per_resouce_block = 12;
%% Transmitter Parameters%%
TX_Para.qam_modulation_order = 16; % 4,16,64,128
TX_Para.ofdm_symbol_duration_ms = (1 / SYS_Para.subcarrier_spacing_khz);
TX_Para.cp_duration_ms = TX_Para.ofdm_symbol_duration_ms * SYS_Para.cyclic_prefix_ratio;
TX_Para.ofdm_with_cp_duration_ms = TX_Para.ofdm_symbol_duration_ms * (1 + SYS_Para.cyclic_prefix_ratio);
TX_Para.data_subcarrier_num = get_data_subcarrier_num(SYS_Para);
TX_Para.fft_size = get_fft_size(SYS_Para,TX_Para);
TX_Para.cyclic_prefix_length = ceil(TX_Para.fft_size * SYS_Para.cyclic_prefix_ratio);
%% Channel Parameters%%
CH_Para.channel_type = 'TDL_A'; % 'AWGN', 'TDL_A', 'TDL_B', 'TDL_C', 'TDL_D', 'TDL_E'
CH_Para.delay_spread_ns = set_delay_spread_ns(CH_Para,30);
CH_Para.max_doppler_shift_hz = set_max_doppler_shift_hz(CH_Para,500); %300 or inf
CH_Para.channel_fft_size = 2 * TX_Para.fft_size;
CH_Para.sim_time_grid_ms = TX_Para.ofdm_with_cp_duration_ms;         % time-domain unit duration
CH_Para.sim_delay_grid_ns = get_sim_delay_grid_ns(SYS_Para,CH_Para); % delay-domain unit duration
CH_Para.sim_frequency_grid_khz = (1e6 / CH_Para.sim_delay_grid_ns) / CH_Para.channel_fft_size;
%% Receiver Parameters%%
if (CH_Para.max_doppler_shift_hz >= (1e3 * SYS_Para.subcarrier_spacing_khz))
   warning('The setting of subcarrier_spacing_khz may be wrong!')
end
%% Derivation Functions
function data_subcarrier_num = get_data_subcarrier_num(SYS_Para)
    scs_idx = find(SYS_Para.numerology.Bandwidth.SCS_kHz == SYS_Para.subcarrier_spacing_khz);
    read_rb_num = strcat(['SYS_Para.numerology.Bandwidth.BW_' num2str(SYS_Para.signal_bandwidth_mhz) 'MHz(scs_idx)']);
    rb_num = eval(read_rb_num);
    if isnan(rb_num)
       error('The combination (SCS,BW) is wrong!')
    end
    data_subcarrier_num = SYS_Para.subcarriers_per_resouce_block * rb_num;
end

function fft_size = get_fft_size(SYS_Para,TX_Para)
    fft_size = 2^(ceil(log2(TX_Para.data_subcarrier_num)));    
    if (fft_size > SYS_Para.max_fft_size)
       error('The modulation FFT size is unexpected!')
    end
end

function sim_delay_grid_ns = get_sim_delay_grid_ns(SYS_Para,CH_Para)
    sim_delay_grid_ns = 1e6 / (CH_Para.channel_fft_size * (SYS_Para.subcarrier_spacing_khz/2));
end

function [delay_spread_ns] = set_delay_spread_ns(CH_Para,delay_spread_ns_2_set)
    if (strcmp(CH_Para.channel_type,'AWGN'))
        delay_spread_ns = 0;
    else
        delay_spread_ns = delay_spread_ns_2_set;
    end
end

function [max_doppler_shift_hz] = set_max_doppler_shift_hz(CH_Para,max_doppler_shift_hz_2_set)
    if (strcmp(CH_Para.channel_type,'AWGN'))
        max_doppler_shift_hz = 0;
    else
        max_doppler_shift_hz = max_doppler_shift_hz_2_set;
    end
end
