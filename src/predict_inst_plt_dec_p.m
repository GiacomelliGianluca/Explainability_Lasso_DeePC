function fig_inst_dec = predict_inst_plt_dec_p(inst, dec, param, sub_num, sub_ind_H, str_x, str_y, str_varname, plt_opt, fig_index)
% predict_inst_dec plots a decision variable dec for the instants inst
% obtained with the parameter param. 
% dec has the columns as the different instants

% - inst: instants where the dec is plotted
% - dec: decision variable
% - param: parameter
% - sub_num: number of submatrices
% - sub_ind_H: matrix containing the indexes of subHankels (if any)
% - plt_opt: plot options



% Plotting weights in the different time instants
if fig_index ~=0
   fig_inst_dec = figure(fig_index);
   hold on, grid on, 
   col = copper(length(inst));
   plt_w_inst = NaN*zeros(length(inst),1);
   for i = 1 : length(inst)
      plt_w_inst(i) = plot(1:1:length(dec(:,1)), dec(:,inst(i)+1), 'Color', col(i,:), 'linewidth', 1, 'DisplayName', strcat('$t=', num2str(inst(i)),'$'));
   end
   xlim([1 931])
   x_title = xlabel(str_x,'Interpreter','latex');
   y_title = ylabel(str_y,'Interpreter','latex');
   h_title = title(str_varname, 'Interpreter','latex', 'FontWeight', 'bold');
   if sub_num > 0 
      xline(sub_ind_H(:,1), 'Color', [0 0 0], 'LineStyle', '-', 'DisplayName', ' ');
      xline(sub_ind_H(:,2), 'Color', [0 0 0], 'LineStyle', '-', 'DisplayName', ' ');
      plt_south_coor = min(ylim) + (max(ylim) - min(ylim))*0.05;
      for i=1 : sub_num 
          if sub_ind_H(i,3) == 10
              text_str = '1';
          else
              text_str = '2';
          end
          text(sub_ind_H(i,1) + (sub_ind_H(i,2) - sub_ind_H(i,1))/2,plt_south_coor,['OP ',text_str],'FontSize',plt_opt.fontsize, 'FontName', plt_opt.fontname, 'HorizontalAlignment','center')
      end
    end
   legend(plt_w_inst)
   legend('Location','best', 'Interpreter','latex');
   % set(h_title, 'FontSize', plt_opt.title_font_dim);
   set(h_title, 'FontSize', 10);
   set(y_title, 'FontSize', plt_opt.y_font_dim);
   set(x_title, 'FontSize', plt_opt.x_font_dim);
   fig_inst_dec.Units               = plt_opt.plot_unit;
   fig_inst_dec.Position(3)         = plt_opt.plot_dim_x;
   fig_inst_dec.Position(4)         = plt_opt.plot_dim_y;
   set(fig_inst_dec.Children, ...
        'FontName',     plt_opt.fontname, ...
        'FontSize',     plt_opt.fontsize - 2);
   set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
   fig_inst_dec.PaperPositionMode   = 'auto';
end



end

