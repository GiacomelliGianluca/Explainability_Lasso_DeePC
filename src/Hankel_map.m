function fig_H_map = Hankel_map(value_OPs, H_label, plt_opt, fig_index, string_matr_type)
% Function that plots the data stored in the Hankel matrix H or Mosaic matrix
% depending on their value. 

% - H: Hankel matrix
% - M: number of OPs
% - conds: conditions on which we define the PWA regions (struct, where .y: output)
% - H_label: Hankel matrix of the same size of H, containing the labels
% - plt_opt: plot options

% Sorting the OPs  
OPs = sort(value_OPs);


% Allocating the figure (columns on x-axis, raws on y-axis)
fig_H_map = figure(fig_index);
grid on, hold on,
str = 'Column index';
x_title = xlabel(str,'Interpreter','latex');
str_y = 'Raw index';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title(['Data Structure - ', string_matr_type, ' matrix'],'Interpreter','latex', 'FontWeight', 'bold');
xlim([1 931]);
xticks([1 200 400 600 800]);
ylim([1 70]);
yticks([1 20 40 60])
ax = gca;
ax.YDir = 'reverse';
ax.XAxisLocation = 'top';
set(h_title, 'FontSize', plt_opt.title_font_dim);
set(y_title, 'FontSize', plt_opt.y_font_dim);
set(x_title, 'FontSize', plt_opt.x_font_dim);
fig_H_map.Units               = plt_opt.plot_unit;
fig_H_map.Position(3)         = plt_opt.plot_dim_x;
fig_H_map.Position(4)         = plt_opt.plot_dim_y;
set(fig_H_map.Children, ...
    'FontName',     plt_opt.fontname, ...
    'FontSize',     plt_opt.fontsize-2);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_H_map.PaperPositionMode   = 'auto';

% Specify the size
sz = 4;

% Identifying the coordinates in the matrix (you cannot use H, becuase it
% has zeroes and so find does not recognize these entries)
[raws, col] = find(H_label);

% Extraction of the labels
label = H_label(sub2ind(size(H_label), raws, col));


scatter(col(label == 1), raws(label == 1), sz, [0.10196078431372549 0.5019607843137255 0.7333333333333333],'filled'); % data OP_1
scatter(col(label == 2), raws(label == 2), sz, [0.6274509803921569 0 0],'filled');  % data OP_2

% Legenda per le condizioni if-else
legend('OP 1', 'OP 2', 'Location','best');


end

