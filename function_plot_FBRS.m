
function function_plot_FBRS(T, all_fre, FBRS)

figure
imagesc(T, 1:length(all_fre), FBRS);
set(gca, 'YDir', 'normal')
colormap('jet');
colorbar

xlabel('Time/s')
ylabel('Frequency/Hz')

n_ticks = 5; 

y_indices = round(linspace(1, length(all_fre), n_ticks));
yticks(y_indices); 

yticklabels(arrayfun(@(x) sprintf('%.2f', x), all_fre(y_indices), 'UniformOutput', false));

xticks(linspace(min(T), max(T), n_ticks));
xticklabels(arrayfun(@(x) sprintf('%.2f', x), xticks, 'UniformOutput', false));

set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 50)

end