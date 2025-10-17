clc
clear
close all

%%
load chirp

data = y;
fs = Fs;

%%
wavename = 'db1';
frame_length = 256;
frame_stride = 32;
level = 7; 
E_ratio = 0.3; 

[log_S, F] = FBRS(data, fs, frame_length, frame_stride, level, wavename, E_ratio);

%%
num_frames = size(log_S,2);
start_indices = (0:num_frames-1) * frame_stride;
center_indices = start_indices + (frame_length/2);

T = center_indices / fs;

function_plot_FBRS(T, F, log_S)