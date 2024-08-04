classdef ChannelClass < handle
    properties (SetAccess=private)
        SYS_Para
        TX_Para
        CH_Para
                
        overlap_save_para % L:circular convolution length & (P-1): number of overlapping poionts
        
        channel_impulse_response
        channel_grid_ms
        
        channel_profile % Tap_Index, Normalized_Delay, Power_dB, Fading_Distribution
        fader_phases % theta[tap_number], phi[tap_number][kSineWaveNumber], varphi[tap_number][kSineWaveNumber]
    end

    % Parameters
    properties (Dependent)
        tap_number
        coherence_time_ms
        coherence_time_grid
        update_interval_grid % next time to update channel impulse response
        max_delay_ns
        
        delay_grid_index  % scaled tap delay --> delay domain grid index (1 ~ 4096 (channel_fft_size))
        
        power_scaling_factors % normalize channel power to be 1
    end
    
    properties (Constant)
        channel_filename = 'channel_profiles.xlsx';
        kSineWaveNumber = 32; % generation of Rayleigth fading waveforms
        kChannelUpdateFreq = 100; % update_interval  =  coherence_time / kChannelUpdateFreq;
        kEnableAWGN = true;
    end
    
    methods
        function obj = ChannelClass(SYS_Para,TX_Para,CH_Para)
            obj.SYS_Para = SYS_Para;
            obj.TX_Para = TX_Para;
            obj.CH_Para = CH_Para;
            if (strcmp(obj.CH_Para.channel_type,'AWGN'))
                obj.channel_profile = [];
            else
                obj.channel_profile = readtable(obj.channel_filename,'Sheet',obj.CH_Para.channel_type);
            end
            obj.initalize_rayleigh_fader_phases;
            obj.initalize_overlap_save_para;
        end
        
        function initalize_rayleigh_fader_phases(obj)
            rng(0,'twister');
            obj.fader_phases.theta = 2 * pi * rand(obj.tap_number,1) - pi;
            obj.fader_phases.phi = 2 * pi * rand(obj.tap_number,obj.kSineWaveNumber) - pi;
            obj.fader_phases.varphi = 2 * pi * rand(obj.tap_number,obj.kSineWaveNumber) - pi;
        end
        
        function initalize_overlap_save_para(obj)
            obj.overlap_save_para.L = 2 * obj.TX_Para.fft_size;
            obj.overlap_save_para.P = (obj.overlap_save_para.L - (obj.TX_Para.fft_size + obj.TX_Para.cyclic_prefix_length)) + 1;
            obj.overlap_save_para.lastpoints = zeros((obj.overlap_save_para.P-1),1);
        end
        % Dependent Parameters
        function val = get.tap_number(obj)
            if (strcmp(obj.CH_Para.channel_type,'AWGN'))
                val = 1;
            else
                val = length(obj.channel_profile.Tap_Index);
            end
        end
        
        function val = get.coherence_time_ms(obj)
            val = 1e3 * (0.4 / obj.CH_Para.max_doppler_shift_hz);
        end
        
        function val = get.coherence_time_grid(obj)
            % Make sure that the min. coherence_time_grid = 1
            val = max(floor(obj.coherence_time_ms / obj.TX_Para.ofdm_with_cp_duration_ms),1);
        end
        
        function val = get.update_interval_grid(obj)
            % Make sure that the min. update_interval_grid = 1
            val = max(2^floor(log2(obj.coherence_time_grid / obj.kChannelUpdateFreq)),1);
        end
        
        function val = get.max_delay_ns(obj)
            val = obj.CH_Para.delay_spread_ns * obj.channel_profile.Normalized_Delay(end);
        end
        
        function val = get.power_scaling_factors(obj)
            if (strcmp(obj.CH_Para.channel_type,'AWGN'))
                val = 1;
            else
                channel_power = sum(10.^(obj.channel_profile.Power_dB / 10));
                val = sqrt((1/channel_power) * 10.^(obj.channel_profile.Power_dB / 10));
            end
        end
        
        function val = get.delay_grid_index(obj)
            if any(strcmp(obj.CH_Para.channel_type,{'AWGN','Single_Path_Rayleigh'}))
                val = 1;
            else
                val = floor((obj.CH_Para.delay_spread_ns * obj.channel_profile.Normalized_Delay) ./ obj.CH_Para.sim_delay_grid_ns) + 1;
            end
        end
        
        % Main Functions
        function [cir_tap_time,update_interval_grid] = generate_cir_tap_time(obj,time_interval_to_sim)
            time_length = time_interval_to_sim(end) - time_interval_to_sim(1) + 1;
            cir_tap_time = zeros(obj.tap_number,time_length);
            if (strcmp(obj.CH_Para.channel_type,'AWGN'))
                cir_tap_time = ones(1,time_length);
            elseif ismember(obj.CH_Para.channel_type,{'TDL_A','TDL_B','TDL_C'})  
                if (obj.CH_Para.max_doppler_shift_hz == inf)
                    cir_tap_time = sqrt(1/2) * (randn(obj.tap_number,time_length) + 1i * randn(obj.tap_number,time_length));
                else
                    update_time_ms = obj.TX_Para.ofdm_with_cp_duration_ms * (ceil((time_interval_to_sim(1):time_interval_to_sim(end))/obj.update_interval_grid) - 1);
                    for tap_idx = 1:obj.tap_number
                        for wave_idx = 1:obj.kSineWaveNumber
                            cir_tap_time(tap_idx,:) = cir_tap_time(tap_idx,:) ...
                                                    + cos(2*pi*obj.CH_Para.max_doppler_shift_hz*cos((2*pi*wave_idx-pi+obj.fader_phases.theta(tap_idx))/(4*obj.kSineWaveNumber))*1e-3*update_time_ms + obj.fader_phases.phi(tap_idx,wave_idx)) ...
                                                 + 1i*cos(2*pi*obj.CH_Para.max_doppler_shift_hz*sin((2*pi*wave_idx-pi+obj.fader_phases.theta(tap_idx))/(4*obj.kSineWaveNumber))*1e-3*update_time_ms + obj.fader_phases.varphi(tap_idx,wave_idx));
                        end
                    end
                    cir_tap_time = sqrt(2/obj.kSineWaveNumber) * cir_tap_time;
                end
            elseif ismember(obj.CH_Para.channel_type,{'TDL_D','TDL_E'})
            else
                error('The channel_type is wrong!')
            end
            update_interval_grid = obj.update_interval_grid;
        end
        
        function [cir_delay_time,delay_grid_ns,update_interval_grid] = generate_power_scaled_cir_delay_time(obj,time_interval_to_sim)
            [cir_tap_time,update_interval_grid] = obj.generate_cir_tap_time(time_interval_to_sim);
            cir_tap_time_power_scaled = diag(obj.power_scaling_factors) * cir_tap_time;
%             normalized_channel_power = sum(sum(abs(real(cir_tap_time_power_scaled)).^2,1))/length(cir_tap_time);
            cir_delay_time = zeros(obj.CH_Para.channel_fft_size,size(cir_tap_time,2));
            % Map tap to delay according to its index
            for tap_idx = 1:length(obj.delay_grid_index)
                cir_delay_time(obj.delay_grid_index(tap_idx),:) = cir_delay_time(obj.delay_grid_index(tap_idx),:) + cir_tap_time_power_scaled(tap_idx,:);
            end
            delay_grid_ns = obj.CH_Para.sim_delay_grid_ns;      
        end
        
        function [cir_freq_time,freq_grid_khz,update_interval_grid] = generate_cir_freq_time(obj,time_interval_to_sim)
            [cir_delay_time,~,update_interval_grid] = obj.generate_power_scaled_cir_delay_time(time_interval_to_sim);
            cir_freq_time = fft(cir_delay_time);
            freq_grid_khz = obj.CH_Para.sim_frequency_grid_khz; 
        end
        
        function [low_pass_channel_output,cir_freq_time] = get_low_pass_channel_output(obj,snr_db,modulated_low_pass_signal,time_interval_to_sim)
            [cir_freq_time,~,~] = obj.generate_cir_freq_time(time_interval_to_sim);
            % overlap save method
            signal_mat = reshape(modulated_low_pass_signal,(obj.TX_Para.fft_size + obj.TX_Para.cyclic_prefix_length),[]);
            signal_mat_overlap_save = zeros(obj.overlap_save_para.L,size(signal_mat,2));
            for column_idx = 1:size(signal_mat,2)
                if (column_idx == 1)
                    signal_mat_overlap_save(:,column_idx) = [obj.overlap_save_para.lastpoints;signal_mat(:,column_idx)];
                else
                    signal_mat_overlap_save(:,column_idx) = [signal_mat((end-obj.overlap_save_para.P+2):end,column_idx-1);signal_mat(:,column_idx)];
                end
            end
            
            signal_mat_overlap_save_f = fft(signal_mat_overlap_save);
            
            ch_out_mat_f = cir_freq_time(1:obj.overlap_save_para.L,:) .* signal_mat_overlap_save_f;
            
            ch_out_mat_t = ifft(ch_out_mat_f);
            
            ch_out_mat_t = ch_out_mat_t(obj.overlap_save_para.P:end,:);
            
            % 2D --> 1D sequence
            ch_out_t = ch_out_mat_t(:);
            
            % Additive Noise
            % Assume that conv(sig,ch) sequence has normalized power = 1.
            noise_pwr = 10^(-1 * snr_db / 10); 
            awgn = sqrt(noise_pwr/2) * (randn(length(ch_out_t),1) + 1i * randn(length(ch_out_t),1));

            if (obj.kEnableAWGN)
                low_pass_channel_output = ch_out_t + awgn;
            else
                low_pass_channel_output = ch_out_t;
            end
        end
        
        % Check Results
        function [channel_power,delay_mean,delay_spread] = check_tdl_channel(obj)
            channel_power = sum(10.^(obj.channel_profile.Power_dB / 10));
            delay_mean = sum(obj.channel_profile.Normalized_Delay .* 10.^(obj.channel_profile.Power_dB / 10)) / channel_power;
            delay_spread = sqrt(sum((obj.channel_profile.Normalized_Delay - delay_mean).^2 .* 10.^(obj.channel_profile.Power_dB / 10)) / channel_power);
        end
        
        function check_correlation_results(obj)
            [cir_tap_time,update_interval] = obj.generate_cir_tap_time([1 2^13]);
            figure
            for tap_idx = 1:size(cir_tap_time,1)
                [acf,lags] = autocorr(real(cir_tap_time(tap_idx,:)),floor(size(cir_tap_time,2)/2));
                time_index = obj.CH_Para.max_doppler_shift_hz * 1e-3 * obj.TX_Para.ofdm_with_cp_duration_ms * update_interval * lags;
                plot(time_index,acf);hold on;
            end
            xlim([0 10])
            xlabel('$f_D*\tau$','Interpreter','latex');ylabel('Autocorrelation of In-Phase Components')
            set(gca,'FontSize',14);set(gcf,'color','w');grid on;
            
            figure
            for tap_idx = 1:size(cir_tap_time,1)
                [acf,lags] = autocorr(imag(cir_tap_time(tap_idx,:)),NumLags=floor(size(cir_tap_time,2)/2));
                time_index = obj.CH_Para.max_doppler_shift_hz * 1e-3 * obj.TX_Para.ofdm_with_cp_duration_ms * update_interval * lags;
                plot(time_index,acf);hold on;
            end
            xlim([0 10])
            xlabel('$f_D*\tau$','Interpreter','latex');ylabel('Autocorrelation of Quadrature Components')
            set(gca,'FontSize',14);set(gcf,'color','w');grid on;
        end
    end
end