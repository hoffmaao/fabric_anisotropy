function h_fig = get_figures(num_figs, visible_flag, user_data)
% Stub of OPR get_figures for standalone testing (always invisible)
h_fig = [];
for k = 1:num_figs
  h_fig(end+1) = figure('Visible','off');
end
end
