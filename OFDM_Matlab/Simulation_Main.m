clear;clc;close all;

%% Create an OFDM System
% Load Parameters 
run Parameters.m
%% AWGN
OFDM_TX = TransmitterClass(SYS_Para,TX_Para);
CH = ChannelClass(SYS_Para,TX_Para,CH_Para);
CH.check_correlation_results();
OFDM_RX = ReceiverClass(SYS_Para,TX_Para,CH_Para);
snr_db_range = 0:2:20;
Result = StatisticClass();
Result.plot_ser_vs_snr(OFDM_TX,CH,OFDM_RX,snr_db_range);
