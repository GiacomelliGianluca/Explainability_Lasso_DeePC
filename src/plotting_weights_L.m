function fig_weights= plotting_weights_L(weights, sub_num, sub_ind_H, plt_opt, fig_index)
% plotting_weights_L plots the different weights selected by DeePC over the
% prediction horizon

% - weights: matrix of weights, where a column represent the weights at a
%            certain instant
% - sub_num: number of submatrices
% - sub_ind_H: matrix containing the indexes of subHankels (if any)
% - plt_opt: plot options

% - fig_weights: plot of the weights in the different time instants

str_x = '$k$ [column]';
str_y = '$[g]_k$';

% Plotting weights in the different time instants
if fig_index ~=0
   fig_weights = figure(fig_index);
   hold on, grid on, 
   col = copper(length(weights(1,:)));
   % plt_w: weight plots
   plt_w = NaN*zeros(length(weights(1,:)),1);
   for i = 1 : length(weights(1,:))
      plt_w(i) = plot(1:1:length(weights(:,1)), weights(:,i), 'Color', col(i,:), 'DisplayName', ['Instant ', num2str(i)]);
   end
   x_title_w1 = xlabel(str_x,'Interpreter','latex');
   y_title_w1 = ylabel(str_y,'Interpreter','latex');
   str_title = 'Data Selection';
   h_title_w1 = title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');
   
   if sub_num > 0 
      xline(sub_ind_H(:,1), 'Color', [0 0 0], 'LineStyle', '--', 'DisplayName', ' ');
      xline(sub_ind_H(:,2), 'Color', [0 0 0], 'LineStyle', '--', 'DisplayName', ' ');
      plt_south_coor = min(ylim) + (max(ylim) - min(ylim))*0.05;
      for i=1 : sub_num 
          text(sub_ind_H(i,1) + (sub_ind_H(i,2) - sub_ind_H(i,1))/2,plt_south_coor,['OP ', num2str(sub_ind_H(i,3))],'FontSize',plt_opt.fontsize, 'FontName', plt_opt.fontname, 'HorizontalAlignment','center')
      end
   end
   legend(plt_w)
   legend( 'Location','eastoutside');
   set(h_title_w1, 'FontSize', plt_opt.title_font_dim);
   set(y_title_w1, 'FontSize', plt_opt.y_font_dim);
   set(x_title_w1, 'FontSize', plt_opt.x_font_dim);
   fig_weights.Units               = plt_opt.plot_unit;
   fig_weights.Position(3)         = plt_opt.plot_dim_x;
   fig_weights.Position(4)         = plt_opt.plot_dim_y;
   set(fig_weights.Children, ...
        'FontName',     plt_opt.fontname, ...
        'FontSize',     plt_opt.fontsize - 2);
   set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
   fig_weights.PaperPositionMode   = 'auto';
end


end