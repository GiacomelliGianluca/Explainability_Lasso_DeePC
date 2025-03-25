function [sub_ind_data, sub_ind_H, sub_num, fig_H_sub] = Hankel_sub(H_y, H_u, value_OPs, N_h, n_x, n_u, conds, fig_H, plt_opt, fig_index)
% Function that identifies the submatrixes associated to the different
% OPs in the Hankel matrix H 

% - sub_ind_data: matrix of the indexes of the submatrices associated with the data
% (where the data considered are the ones of the first column and the
% last raw (!! hence, with mosaic not all data are covered and there is not
% an alignment with the database indexes!!)). See below for detail
% - sub_ind_H: matrix of the indexes associated with the position (raw and column)
% in the Hankel matrix of the submatrices. See below for detail
% - sub_num: number of subamtrices
% - fig_H_sub: handle for the submatrix figure

% - H: Hankel matrix
% - value_OPs: column vector with the values of the OPs (in deg)
% - N_h: Hankel matrix depth
% - n_x: state degree
% - n_u: input degree
% - conds: conditions on which we define the partitions (struct, where .y: output)
% - fig_H: Hankel (or Mosaic) figure handle

% Sorting the OPs  
OPs = sort(value_OPs);
% n_p: number of OPs
n_p = length(OPs);

% Subspace predictor id conditions
SP_id_cond = (N_h + n_x)*(n_u + 1) - 1;

% Extracting the switching conditions [rad]
cond_y = pi/180*conds.y;


% Retrieving the data vector
data_y = [H_y(:,1)'   H_y(end,2:end)];
data_u = [H_u(:,1)'   H_u(end,2:end)];

% Indexes initialization
% -Column 1: starting index -Column 2: ending index -Column 3: belonging OP
% Each raw describe a differen submatrices
sub_ind_data = [];

% f_b_submatr: flag identifying if a submatrix is being identify
f_b_submatr = 0;

% sub_num: counter for the number of submatrices
sub_num = 0;

for i = 1 : (length(data_y) - SP_id_cond + 1 )
    
    % Identification of the submatrices associated to the OP 1
    if (sum(data_y(i : (i - 1) + SP_id_cond) <= cond_y) == SP_id_cond) || (f_b_submatr == 1)
       if f_b_submatr == 0
           sub_num = sub_num + 1;
           f_b_submatr = 1;
           sub_ind_data(sub_num, 1) = i;
           sub_ind_data(sub_num, 3) = OPs(1);
       end
       %
       if sum(data_y(i : (i - 1) + SP_id_cond) <= cond_y) ~= SP_id_cond
           f_b_submatr = 0;
           sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
       end
       %
       if i == (length(data_y) - SP_id_cond + 1) && (f_b_submatr == 1)
           f_b_submatr = 0;
           sub_ind_data(sub_num, 2) = length(data_y);
       end
    end

    % Identification of the submatrices associated to the OP n_m
    if (sum(data_y(i : (i - 1) + SP_id_cond) > cond_y) == SP_id_cond) || (f_b_submatr == n_p)
       if f_b_submatr == 0
           sub_num = sub_num + 1;
           f_b_submatr = n_p;
           sub_ind_data(sub_num, 1) = i;
           sub_ind_data(sub_num, 3) = OPs(end);
       end
       %
       if sum(data_y(i : (i - 1) + SP_id_cond) > cond_y) ~= SP_id_cond
           f_b_submatr = 0;
           sub_ind_data(sub_num, 2) = (i - 2) + SP_id_cond;
       end
    end


   % Exit case
   if i == (length(data_y) - SP_id_cond + 1) && (f_b_submatr > 0)
       f_b_submatr = 0;
       sub_ind_data(sub_num, 2) = length(data_y);
   end

end

% Checking persistency of excitation of the submatrices
for j = 1 : sub_num
   PE_check(data_u(sub_ind_data(j,1):sub_ind_data(j,2)), (n_u* (N_h + n_x)));
end

% Plotting (the final indexes)
if sub_num > 0 
    if fig_index~=0
        % Allocating the figure (columns on x-axis, raws on y-axis)
        fig_H_sub = figure(fig_index); hold on
        CloneFig(fig_H, fig_H_sub);
        ax = gca;
        ax.YDir = 'reverse';
        ax.XAxisLocation = 'top';
        fig_H_sub.Units               = plt_opt.plot_unit;
        fig_H_sub.Position(3)         = plt_opt.plot_dim_x;
        fig_H_sub.Position(4)         = plt_opt.plot_dim_y;
        set(fig_H_sub.Children, ...
            'FontName',     plt_opt.fontname, ...
            'FontSize',     plt_opt.fontsize-2);
        set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
        fig_H_sub.PaperPositionMode   = 'auto';
        % Adjusting the legend
        allChildren = get(gca, 'Children');                % list of all objects on axes
        displayNames = get(allChildren, 'DisplayName');    % list of all legend display names
        delete(allChildren(1:length(displayNames)/2));
        grai = gray(sub_num);
        for j = 1 : sub_num
           xline(sub_ind_data(j,1), 'DisplayName', ['Submatrix ', num2str(j),' - begin'], 'Color', grai(j,:), LineWidth=2)
           xline(sub_ind_data(j,2) - N_h + 1, 'DisplayName', ['Submatrix ', num2str(j),' - end'], 'Color', grai(j,:), LineWidth=2)
        end
        legend('NumColumns',(2*n_p)-1)
    else
        fig_H_sub = NaN;
    end

    % Identify Hankel indexes 
    sub_ind_H(:,1) = sub_ind_data(:, 1);
    sub_ind_H(:,3) = sub_ind_data(:, 3);
    for k = 1 : sub_num
        sub_ind_H(k, 2) = sub_ind_data(k, 2) - N_h + 1;
    end
else 
    fig_H_sub = NaN;
    sub_ind_H = sub_ind_data;
end

end

