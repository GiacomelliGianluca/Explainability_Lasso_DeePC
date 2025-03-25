function fig_handle = dim3plotting(varX, varY, varZ, title_string, x_string, y_string, z_string, plt_opt, fig_index)
% dim3plotting takes the variables VarX, VarY, and VarZ and plots them.
% It is assumed that there is a positional correspondance in the vectors of
% the variables and the one of varColor

% fig_index: index of the current figure
% varColor_string: prefix string used in the legend for VarColor
% plt_opt: plot options


% fig_handle: handle of the figure


% Plotting the figure
fig_handle = figure(fig_index);
hold on, grid on
x_title = xlabel(x_string,'Interpreter','latex');
y_title = ylabel(y_string,'Interpreter','latex');
z_title = zlabel(z_string,'Interpreter','latex'); 
h_title = title(title_string, 'Interpreter','latex', 'FontWeight', 'bold');

plot3(varX, varY, varZ, 'linewidth',plt_opt.thick, 'Color', '#9e9ac8');
set(gca, 'XScale', 'log')
view(8,15)


set(h_title, 'FontSize', plt_opt.fontsize);
set(z_title, 'FontSize', plt_opt.y_font_dim);
set(y_title, 'FontSize', plt_opt.x_font_dim);
set(x_title, 'FontSize', plt_opt.x_font_dim);
fig_handle.Units               = plt_opt.plot_unit;
fig_handle.Position(3)         = plt_opt.plot_dim_x;
fig_handle.Position(4)         = plt_opt.plot_dim_y;
set(fig_handle.Children, ...
    'FontName',     plt_opt.fontname, ...
    'FontSize',     plt_opt.fontsize - 2);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_handle.PaperPositionMode   = 'auto';
set(get(gca,'ylabel'),'rotation',60)

end

