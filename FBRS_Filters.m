
function [filter_matrix, core_fre, freq_values] = FBRS_Filters(data, wavename, level, E_ratio, fs)

%----------------------------------------------%
% Author：Rui Qin
% paper：Frequency Band Recalibration Spectrogram
%----------------------------------------------%

[~,~,~,~,E] = modwpt(data, level, wavename);

nodes = linspace(2^level-1, 2^(level+1)-2, 2^level);

[sort_E, sort_index] = sort(E, 'descend');

num_components = round(length(sort_E) * E_ratio);
selected_indices = sort_index(1:num_components);

filtered_index_0 = adjust_indices(nodes, selected_indices);

nodes_zeros = Reset_nodes(nodes, filtered_index_0);

node_index = [];
filtered_index = [];
all_index = [];

for k = 1:level-1
    if k==1  
        [new_index_1, filtered_index_1] = Stacking_energy_upwards(nodes_zeros, E, E_ratio);
        node_index = new_index_1;
        filtered_index = filtered_index_1;
        all_index = [all_index filtered_index];
    else   
        nodes_zeros = Reset_nodes(node_index, filtered_index);
        [new_index_2, filtered_index_2] = Stacking_energy_upwards(nodes_zeros, E, E_ratio);
        node_index = new_index_2;
        filtered_index = filtered_index_2;
        all_index = [all_index filtered_index];
    end
end

all_index = [filtered_index_0 all_index];

all_index = sort(all_index);

for i = 1:length(all_index)
    node_level(i,:) = find_level(all_index(:,i), level);
end

fre_index_start = [];
fre_index_end = [];
for i = 1:length(all_index)
    n = node_level(i,:);
    index = all_index(:,i) - (2^n-2);
    fre_res = Filter_frequency_resolution(fs, n);

    fre_start = fre_res(index);
    fre_end = fre_res(index+1);

    fre_index_start = [fre_index_start fre_start];
    fre_index_end = [fre_index_end fre_end];
end

fre_index = [fre_index_start fre_index_end];
fre_index = sort(fre_index);

freq_boundaries = unique(fre_index);

overlap_ratio = 0;
new_freq_boundaries = create_overlapping_boundaries(freq_boundaries, overlap_ratio);

nfft = 1024;
[filter_matrix, core_fre, freq_values, ~] = generate_filter_matrix(new_freq_boundaries, nfft);

end


%%
function nodes_zeros = Reset_nodes(original_nodes, filter_nodes)

zero_idx = ismember(original_nodes, filter_nodes);

original_nodes(zero_idx) = 0;

nodes_zeros = original_nodes;

end


%% 
function [filter_matrix, core_fre, freq_values, num_filters] = generate_filter_matrix(new_freq_boundaries, nfft)

freq_values = linspace(min(new_freq_boundaries), max(new_freq_boundaries), nfft);

num_filters = length(new_freq_boundaries) / 2;
filter_matrix = zeros(num_filters, length(freq_values));


core_fre = [];

for i = 1:num_filters
    start_freq = new_freq_boundaries(2 * i - 1);
    end_freq = new_freq_boundaries(2 * i);
    
    peak_freq = (start_freq + end_freq) / 2;
    core_fre = [core_fre peak_freq ];
    x = [start_freq, peak_freq, end_freq];
    y = [0, 1, 0];
    
    filter_response = interp1(x, y, freq_values, 'linear', 0);

    filter_response = filter_response / max(filter_response);
    
    filter_matrix(i, :) = filter_response;
end

end


%%
function new_freq_boundaries = create_overlapping_boundaries(freq_boundaries, overlap_ratio)

new_freq_boundaries = [];
for i = 1:length(freq_boundaries) - 1
    start_freq = freq_boundaries(i);
    end_freq = freq_boundaries(i + 1);
    overlap = (end_freq - start_freq) * overlap_ratio;
    
    if i == 1
        new_freq_boundaries = [new_freq_boundaries, start_freq, end_freq];
    else
        new_freq_boundaries = [new_freq_boundaries, start_freq - overlap, end_freq];
    end
end

end

%% 

function filtered_index = adjust_indices(nodes, selected_indices)
    filtered_index = [];
    for i = 1:length(selected_indices)
        current_index = nodes(selected_indices(i));
        if mod(current_index, 2) == 1
            filtered_index = [filtered_index, current_index, current_index + 1];
        else
            filtered_index = [filtered_index, current_index - 1, current_index];
        end
    end
    filtered_index = unique(filtered_index);
end


%% 
function n = find_level(index, level)

    for n = 1:level
        if index >= (2^n)-1 && index <= (2^(n+1))-2
            return
        end
    end
    
    n = -1;
end

%% 

function fre_res = Filter_frequency_resolution(fs, level)

fre_res = linspace(0, fs/2, 2^level+1);

end

%% 

function [new_index, filtered_index, new_E] = Stacking_energy_upwards(nodes_zeros, E, E_ratio)

new_E = []; 
new_index = []; 

for i = 1:2:length(nodes_zeros)
    if nodes_zeros(i) ~= 0
        if mod(nodes_zeros(i), 2) == 1
            new_E = [new_E, E(i) + E(i+1)];
            ind = ((nodes_zeros(i) + 1)) / 2 - 1;
            new_index = [new_index, ind];
        else
            new_E = [new_E, E(i) + E(i-1)];
            ind = ((nodes_zeros(i))) / 2 - 1;
            new_index = [new_index, ind];
        end
    end

end

isolated_node = [];
for i = 1:length(new_index)
    if mod(new_index(i), 2) == 0
        if ~ismember(new_index(i) - 1, new_index)
            isolated_node = [isolated_node, new_index(i)];
        end
    else 
        if ~ismember(new_index(i) + 1, new_index)
            isolated_node = [isolated_node, new_index(i)];
        end
    end
end

if isequal(new_index, isolated_node)
    disp('All values in new_index are isolated values; exit the function.');
    filtered_index = new_index;
    return; 
end

[sort_E, sort_index] = sort(new_E, 'descend');

N = length(sort_E);
num_components = ceil(N * E_ratio);
% disp(['num_comp:', num2str(num_components)])
selected_indices = sort_index(1:num_components);

filtered_index = adjust_indices(new_index, selected_indices);

filtered_index = [filtered_index isolated_node];

filtered_index = sort(filtered_index);

end

