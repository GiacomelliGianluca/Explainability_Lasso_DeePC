function fig_data_hl = plotting_dataset_hl(time, value_modes, Data_u, Data_y, conds, str_x, str_u, str_y,  plt_set, fig_index)
% Function that plots the data depending on their value against.
% It is assumed that the modes can be reached sequentially, e.g.:
% ... mode_i-1 --> mode_i --> mode_i+1 ...

% - Data_u: input data (in deg)
% - Data_y: output data (in deg)
% - value_modes: column vector with the values of the modes (in deg)
% - conds: conditions on which we define the trajectory blocks (struct,
% where .y: output) (in deg)

% - fig_data_hl: plots of the dataset trajectory in a subplot

% Sorting the modes  
modes = sort(value_modes);

% Specify the size of plot points
sz = 4;
plt_set.plot_dim_x = 8.5;
plt_set.plot_dim_y = 11;

% Allocating the figure (columns on x-axis, raws on y-axis)
fig_data_hl = figure(fig_index);
subplot(2,1,1)
grid on, hold on,
x_title = xlabel(str_x,'Interpreter','latex');
y_title = ylabel(str_u,'Interpreter','latex');
h_title = title('Input Dataset','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_data_hl.Units               = plt_set.plot_unit;
fig_data_hl.Position(3)         = plt_set.plot_dim_x;
fig_data_hl.Position(4)         = plt_set.plot_dim_y;
set(fig_data_hl.Children, ...
    'FontName',     'Times', ...
    'FontSize',     11);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_data_hl.PaperPositionMode   = 'auto';
%
subplot(2,1,2)
grid on, hold on,
x_title = xlabel(str_x,'Interpreter','latex');
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Output dataset','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_data_hl.Units               = plt_set.plot_unit;
fig_data_hl.Position(3)         = plt_set.plot_dim_x;
fig_data_hl.Position(4)         = plt_set.plot_dim_y;
set(fig_data_hl.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_data_hl.PaperPositionMode   = 'auto';

 % Extracting the data collection conditions 
cond_y = conds.y;
    
% Data_m1: Data mode 1
Data_m1 = Data_y < cond_y;
% Data_m2: Data mode 2
Data_m2 = Data_y >= cond_y;

% Without middle points
%input
subplot(2,1,1)
scatter(time(Data_m1), Data_u(Data_m1), sz, [0.10196078431372549 0.5019607843137255 0.7333333333333333],'filled'); % data mode_1
scatter(time(Data_m2), Data_u(Data_m2), sz, [0.6274509803921569 0 0],'filled');  % data mode_2
%output
subplot(2,1,2)
scatter(time(Data_m1), Data_y(Data_m1), sz, [0.10196078431372549 0.5019607843137255 0.7333333333333333],'filled'); % data mode_1
scatter(time(Data_m2), Data_y(Data_m2), sz, [0.6274509803921569 0 0],'filled');  % data mode_2

% Legenda per le condizioni if-else
legend(['Mode ', num2str(modes(1))], ['Mode ', num2str(modes(2))], 'Location','best');

end

%OLD
% % mode_tr: modes with transitions between them.
% % e.g. [mode_1 ; mode_12; mode_2]
% modes_tr = [];
% for i = 1 : (length(modes) - 1)
%     modes_tr = [modes_tr; modes(i); (modes(i) + modes(i+1)) / 2];
% end
% modes_tr = [modes_tr; modes(end)];
% 
% % b_mode: ball of the modes and the transition between them 
% b_mode = repmat(struct('Left', [], 'Right', []), 2*length(modes) - 1, 1);
% % First mode: special case, only right ball
% b_mode(1).Right = (modes_tr(1) + modes_tr(2))/2;
% for j = 2 : (2*length(modes) - 2)
%     % left
%     b_mode(j).Left = b_mode(j-1).Right;
%     % right
%     b_mode(j).Right = (modes_tr(j) + modes_tr(j+1))/2;
% end
% % Last mode: special case, only left ball
% b_mode(end).Left = (modes_tr(end-1) + modes_tr(end))/2;

% % Plotting the data depending on the conditions
% color = hsv(2*length(modes) - 1);
% % Specify the size
% sz = 20;
% % mode_str: current mode; Used for the legend
% mode_str = 1;
% % The number of existing regions is double the number modes minus 1
% for k = 1 : (2*length(modes) - 1)
%     % First mode case 
%     if k == 1
%         subplot(2,1,1)
%         scatter(time(Data_y < b_mode(k).Right), Data_u(Data_y <  b_mode(k).Right), sz, color(k,:), 'filled', 'DisplayName', ['data mode_1 (', num2str( round(modes_tr(1), 2)),'^o)']); % data mode_1
%         subplot(2,1,2)
%         scatter(time(Data_y < b_mode(k).Right), Data_y(Data_y <  b_mode(k).Right), sz, color(k,:), 'filled', 'DisplayName', ['data mode_1 (', num2str( round(modes_tr(1), 2)),'^o)']); % data mode_1
%     else
%         % Last mode case 
%         if k == (2*length(modes) - 1)
%             subplot(2,1,1)
%             scatter(time(Data_y >= b_mode(end).Left), Data_u(Data_y >= b_mode(end).Left), sz, color(k,:), 'filled', 'DisplayName', ['data mode_', num2str(length(modes)), ' (', num2str( round(modes_tr(end), 2)),'^o)']);  % data last mode
%             subplot(2,1,2)
%             scatter(time(Data_y >= b_mode(end).Left), Data_y(Data_y >= b_mode(end).Left), sz, color(k,:), 'filled', 'DisplayName', ['data mode_', num2str(length(modes)), ' (', num2str( round(modes_tr(end), 2)),'^o)']);  % data last mode
%         else
%             % General mode case (if-else condition is for the plotting)
%             if rem(k,2) == 0
%                 % Transition
%                 mode_str = mode_str + 1;   
%                 subplot(2,1,1)
%                 scatter(time((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), Data_u((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), sz, color(k,:), 'filled', 'DisplayName', ['data between mode_', num2str(mode_str - 1), ' and mode_', num2str(mode_str)]); % data mode_1, mode_2
%                 subplot(2,1,2)
%                 scatter(time((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), Data_y((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), sz, color(k,:), 'filled', 'DisplayName', ['data between mode_', num2str(mode_str - 1), ' and mode_', num2str(mode_str)]); % data mode_1, mode_2
%             else
%                 % Mode
%                 subplot(2,1,1)
%                 scatter(time((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), Data_u((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), sz, color(k,:), 'filled', 'DisplayName', ['data mode_', num2str(mode_str), ' (', num2str( round(modes_tr(k), 2)),'^o)']); % data mode_1, mode_2
%                 subplot(2,1,2)
%                 scatter(time((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), Data_y((Data_y >= b_mode(k).Left) & (Data_y < b_mode(k).Right)), sz, color(k,:), 'filled', 'DisplayName', ['data mode_', num2str(mode_str), ' (', num2str( round(modes_tr(k), 2)),'^o)']); % data mode_1, mode_2
%             end
%         end
%     end
% end
% 
% 
% % Legenda per le condizioni if-else
% legend show;
% % legend('Location','eastoutside');

