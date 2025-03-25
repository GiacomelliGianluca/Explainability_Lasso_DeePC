function fig_data = plotting_dataset(time, Data_u, Data_y, str_x, str_u, str_y,  plt_set, fig_index)
% plotting_dataset plots the dataset trajectory Data_u, Data_y in a subplot

% - Data_u: input data
% - Data_y: output data

% - fig_data: plots of the dataset trajectory in a subplot

plt_set.plot_dim_x = 8.5;
plt_set.plot_dim_y = 11;

fig_data = figure(fig_index);
subplot(2,1,1)
plot(time, Data_u,'k','linewidth',plt_set.thick), hold on, grid on, 
x_title = xlabel(str_x,'Interpreter','latex');
y_title = ylabel(str_u,'Interpreter','latex');
h_title = title('Input Dataset','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_data.Units               = plt_set.plot_unit;
fig_data.Position(3)         = plt_set.plot_dim_x;
fig_data.Position(4)         = plt_set.plot_dim_y;
set(fig_data.Children, ...
    'FontName',     'Times', ...
    'FontSize',     11);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_data.PaperPositionMode   = 'auto';
%
subplot(2,1,2)
plot(time, Data_y,'k','linewidth',plt_set.thick), hold on, grid on, 
x_title = xlabel(str_x,'Interpreter','latex');
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Output dataset','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_data.Units               = plt_set.plot_unit;
fig_data.Position(3)         = plt_set.plot_dim_x;
fig_data.Position(4)         = plt_set.plot_dim_y;
set(fig_data.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_data.PaperPositionMode   = 'auto';

end

