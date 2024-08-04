classdef ReceiverClass < handle
    properties (SetAccess=private)
        SYS_Para
        TX_Para
        CH_Para
    end
    
    % Parameters
    properties (Dependent)
    end
    
    properties (Constant)
    end
    
    methods
        function obj = ReceiverClass(SYS_Para,TX_Para,CH_Para)
            obj.SYS_Para = SYS_Para;
            obj.TX_Para = TX_Para;
            obj.CH_Para = CH_Para;
        end
        
        function [data_qam_f_equalized, data_qam_index_detect] = get_data_qam_index_detect(obj,received_low_pass_signal,cir_freq_time)
            % Remove Cyclic Prefix
            data_with_cp_time = reshape(received_low_pass_signal,(obj.TX_Para.fft_size + obj.TX_Para.cyclic_prefix_length),[]);
            data_time = data_with_cp_time((obj.TX_Para.cyclic_prefix_length+1):end,:);
            % FFT Transform
            data_qam_f = fft(data_time) / sqrt(obj.TX_Para.fft_size);
            % Frequency Domain Equalizer
            data_qam_f = data_qam_f(1:obj.TX_Para.data_subcarrier_num,:);
            data_qam_f_equalized = data_qam_f ./ cir_freq_time(1:2:(2*obj.TX_Para.data_subcarrier_num-1),:);
            % QAM Demapping
            data_qam_index_detect = qamdemod(data_qam_f_equalized,obj.TX_Para.qam_modulation_order,'UnitAveragePower',true);
        end
    end
end