function [log_S, core_fre, freq_values] = FBRS(data, fs, frame_length, frame_stride, level, wavename, E_ratio)

[filter_matrix, core_fre, freq_values] = FBRS_Filters(data, wavename, level, E_ratio, fs);

Data_E = FBRS_pre_process(data, frame_length, frame_stride);

S = Data_E * filter_matrix';

log_S = FBRS_log(S);

end
