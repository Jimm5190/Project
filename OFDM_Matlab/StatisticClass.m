classdef StatisticClass < handle
    properties (SetAccess=private)        
        BUFFER_SER_VS_SNR_dB_OVER_DOPPLER_SHIFT
    end
    properties (Dependent)
    end
    properties (Constant)
        kSegmentLengthOfdmSymbol = 1024;
        kNumSegment = 10;
    end
    methods
        function obj = StatisticClass()
        end
        
        function [ser_vs_snr] = evaluate_ser_vs_snr(obj,TX,CH,RX,snr_db_range)
            ser = zeros(length(snr_db_range),1);
            ser_analytical = zeros(length(snr_db_range),1);
            for snr_db = snr_db_range
                symbol_error_count = 0;
                symbol_count = 0;
                for seg_idx = 1:obj.kNumSegment
                    disp(['SNR (dB) = ' num2str(snr_db) ', Segment Number: ' num2str(seg_idx)]);
                    time_interval_to_sim = [((seg_idx-1) * obj.kSegmentLengthOfdmSymbol + 1) (seg_idx * obj.kSegmentLengthOfdmSymbol)];
                    [~, data_qam_index, modulated_low_pass_signal] = TX.get_modulated_low_pass_signal(time_interval_to_sim);
                    [low_pass_channel_output,cir_freq_time] = CH.get_low_pass_channel_output(snr_db,modulated_low_pass_signal,time_interval_to_sim);
                    [~, data_qam_index_detect] = RX.get_data_qam_index_detect(low_pass_channel_output,cir_freq_time);
                    symbol_error_count = symbol_error_count + nnz(data_qam_index - data_qam_index_detect);
                    symbol_count = symbol_count + size(data_qam_index,1) * size(data_qam_index,2); 
                end
                ser(snr_db_range==snr_db) = symbol_error_count / symbol_count;
                ser_analytical(snr_db_range==snr_db) = obj.get_analytical_ser_of_M_QAM(snr_db,TX.TX_Para.qam_modulation_order,CH);
            end
            ser_vs_snr.ser = ser;
            ser_vs_snr.ser_analytical = ser_analytical;
            ser_vs_snr.snr_db = snr_db_range;
        end
        
        function plot_ser_vs_snr(obj,TX,CH,RX,snr_db_range)
            ser_vs_snr = obj.evaluate_ser_vs_snr(TX,CH,RX,snr_db_range);
            figure()
            if (strcmp(CH.CH_Para.channel_type,'AWGN'))
                semilogy(ser_vs_snr.snr_db,ser_vs_snr.ser,'bo','linewidth',2);hold on;
                semilogy(ser_vs_snr.snr_db,ser_vs_snr.ser_analytical,'k--','linewidth',2);grid minor;
                title([num2str(TX.TX_Para.qam_modulation_order) ' QAM, ' CH.CH_Para.channel_type ' Channel'], 'Interpreter', 'none');
                legend('Simulated','Analytical')
            elseif (ismember(CH.CH_Para.channel_type,{'TDL_A','TDL_B','TDL_C'}) && (CH.CH_Para.max_doppler_shift_hz == inf))
                semilogy(ser_vs_snr.snr_db,ser_vs_snr.ser,'bo','linewidth',2);hold on;
                semilogy(ser_vs_snr.snr_db,ser_vs_snr.ser_analytical,'k--','linewidth',2);grid minor;
                title({[num2str(TX.TX_Para.qam_modulation_order) ' QAM, ' CH.CH_Para.channel_type ' Channel, Delay Spread ' num2str(CH.CH_Para.delay_spread_ns) ...
                       ' ns'],'Changing Independently Every OFDM Symbol'}, 'Interpreter', 'none');
                legend('Simulated','Analytical')
            else
                semilogy(ser_vs_snr.snr_db,ser_vs_snr.ser,'b-','linewidth',2);grid minor;
                title([num2str(TX.TX_Para.qam_modulation_order) ' QAM, ' CH.CH_Para.channel_type ' Channel, Delay Spread ' num2str(CH.CH_Para.delay_spread_ns) ...
                       ' ns, Max. Doppler Shift ' num2str(CH.CH_Para.max_doppler_shift_hz) 'Hz'], 'Interpreter', 'none');
                legend('Simulated')
            end
            xlabel('SNR (dB)');ylabel('Symbol Error Rate');
            set(gca,'FontSize', 18);set(gcf,'color','w');
        end      
        
        function [ser_vs_snr_over_doppler_shift] = evaluate_ser_vs_snr_over_doppler_shift(obj,SYS_Para,TX_Para,CH_Para,snr_db_range,max_doppler_shift_hz_range)
            ser_vs_snr_over_doppler_shift.ser = zeros(length(snr_db_range),length(max_doppler_shift_hz_range));
            ser_vs_snr_over_doppler_shift.snr_db = snr_db_range(:);
            for doppler_idx = 1:length(max_doppler_shift_hz_range)
                tic
                CH_Para.max_doppler_shift_hz = max_doppler_shift_hz_range(doppler_idx);
                TX = TransmitterClass(SYS_Para,TX_Para);
                CH = ChannelClass(SYS_Para,TX_Para,CH_Para);
                RX = ReceiverClass(SYS_Para,TX_Para,CH_Para);
                [ser_vs_snr] = obj.evaluate_ser_vs_snr(TX,CH,RX,snr_db_range);
                ser_vs_snr_over_doppler_shift.ser(:,doppler_idx) = ser_vs_snr.ser(:);
                toc
            end
        end
        
        function plot_ser_vs_snr_over_doppler_shift(obj,SYS_Para,TX_Para,CH_Para,snr_db_range,max_doppler_shift_hz_range)
            [ser_vs_snr_over_doppler_shift] = obj.evaluate_ser_vs_snr_over_doppler_shift(SYS_Para,TX_Para,CH_Para,snr_db_range,max_doppler_shift_hz_range);
            figure()
            if (strcmp(CH_Para.channel_type,'AWGN'))
                error('Wrong Channel Type!');
            else
                for doppler_idx = 1:length(max_doppler_shift_hz_range)
                    semilogy(snr_db_range,ser_vs_snr_over_doppler_shift.ser(:,doppler_idx),'--','linewidth',2);hold on;
                    my_legend{doppler_idx} = ['$f_{Max Doppler} =$ ' num2str(max_doppler_shift_hz_range(doppler_idx)) ' Hz'];
                end
                title([num2str(TX_Para.qam_modulation_order) ' QAM, '...
                       'SCS ' num2str(SYS_Para.subcarrier_spacing_khz) 'kHz, '...
                       CH_Para.channel_type ' Channel, Delay Spread ' num2str(CH_Para.delay_spread_ns) ' ns'], 'Interpreter', 'none');
                legend(my_legend, 'Interpreter', 'latex');grid minor;
            end
            xlabel('SNR (dB)');ylabel('Symbol Error Rate');
            set(gca,'FontSize', 18);set(gcf,'color','w');
        end
        
        function [ser] = get_analytical_ser_of_M_QAM(~,snr_db,qam_modulation_order,CH)
            snr_linear = 10^(snr_db/10);
            M = qam_modulation_order;
            if ((strcmp(CH.CH_Para.channel_type,'AWGN')) && (M == 16))
                ser = 3 * qfunc(sqrt(snr_linear/5)) - (9/4) * (qfunc(sqrt(snr_linear/5)))^2;
            elseif (ismember(CH.CH_Para.channel_type,{'TDL_A','TDL_B','TDL_C'}) && (CH.CH_Para.max_doppler_shift_hz == inf))
                g_QAM = 3 / (2 * (M - 1));
                ser = 2 * (1 - 1/sqrt(M)) * (1 - sqrt((g_QAM*snr_linear) / (1 + g_QAM*snr_linear))) + ...
                      (1 - 1/sqrt(M))^2 * ((4/pi) * sqrt((g_QAM*snr_linear) / (1 + g_QAM*snr_linear)) * atan(sqrt((1 + g_QAM*snr_linear) / (g_QAM*snr_linear))) - 1);
            else
                ser = [];
            end
        end
    end
end