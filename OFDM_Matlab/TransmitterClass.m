classdef TransmitterClass < handle
    properties (SetAccess=private)
        SYS_Para
        TX_Para
    end
    % Parameters
    properties (Dependent)        
        sample_duration_ms
        ofdm_symbol_duration_ms
        cyclic_prefix_duration_ms
    end
    
    properties (Constant)

    end
    
    methods
        function obj = TransmitterClass(SYS_Para,TX_Para)
            obj.SYS_Para = SYS_Para;
            obj.TX_Para = TX_Para;
        end
        
        %% Dependent Parameters
        function val = get.sample_duration_ms(obj)
            val = 1 / (1e3 * obj.SYS_Para.signal_bandwidth_mhz);
        end
        
        function val = get.ofdm_symbol_duration_ms(obj)
            val = 1 / obj.SYS_Para.subcarrier_spacing_khz;
        end
        
        function val = get.cyclic_prefix_duration_ms(obj)
            val = (obj.TX_Para.cyclic_prefix_length / obj.TX_Para.fft_size) * obj.ofdm_symbol_duration_ms;
        end
        
        %% Main Functions
        function [data_qam, data_qam_index, modulated_low_pass_signal] = get_modulated_low_pass_signal(obj,time_interval_to_sim )  
            time_interval = time_interval_to_sim(end) - time_interval_to_sim(1) + 1;
            % Generate Random Bits
            data_qam_index = randi([0 (obj.TX_Para.qam_modulation_order-1)],obj.TX_Para.data_subcarrier_num,time_interval); %[fft_size][kOFDMSymbolLength]
           
            % QAM Mapping
            data_qam = qammod(data_qam_index,obj.TX_Para.qam_modulation_order,'UnitAveragePower',true); % [fft_size][kOFDMSymbolLength]
            
            % IFFT
            data_qam_ext = [data_qam; zeros((obj.TX_Para.fft_size-obj.TX_Para.data_subcarrier_num),time_interval)];
            data_time = (conj(dftmtx(obj.TX_Para.fft_size))/sqrt(obj.TX_Para.fft_size))*data_qam_ext;
            
            % Add Cyclic Prefix
            data_with_cp_time = [data_time(end-(obj.TX_Para.cyclic_prefix_length-1):end,:); data_time];
            
            % Parallel to Serial
            modulated_low_pass_signal = data_with_cp_time(:);
            
            % DAC (Omitted)
            
            % Modulated Signal at Carrier Frequency (Omitted)
            
        end
    end
end