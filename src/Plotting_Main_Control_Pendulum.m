% Main script for plotting DDPC reference tracking of a given trajectory
% using Hankel data structure

close all 
clear 
clc

% Before Execute the following, run the script dataset creation 

%%  Plot settings
% plt_set: Plot settings
% Settings
plt_set.y_font_dim = 16;
plt_set.x_font_dim = 11;

plt_set.title_font_dim = 14;

% %thickness
plt_set.thick = 1.6;
% Figure dimensions (inluding frames)
plt_set.plot_dim_x = 12;
plt_set.plot_dim_y = 5;
% Unit of measure
% plt_set.plot_unit = 'centimeters';
plt_set.plot_unit = 'inches'; %used in ppt
% Used font
plt_set.fontname = 'Times';
plt_set.fontsize = 11;

%% Loading Simulation 

% Datasets directory HERE
dataset_dir = [pwd, '\Simulations\'];

% test traj10_const_10_90_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_90_random_traj10_const_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_case2_traj10_const_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_SPC_ok_traj10_const_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_two_blocks_traj10_const_10_90_ini_mode_1.mat';
% 
% test traj10_const_10_170_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_170_random_traj10_const_10_170_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_case2_traj10_const_10_170_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_SPC_ok_traj10_const_10_170_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_two_blocks_traj10_const_10_170_ini_mode_1.mat';
% 
% test traj90_const_90_10_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_90_random_traj90_const_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_case2_traj90_const_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_SPC_ok_traj90_const_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_two_blocks_traj90_const_90_10_ini_mode_1.mat';
%
% test traj170_const_170_10_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_170_random_traj170_const_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_case2_traj170_const_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_SPC_ok_traj170_const_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_two_blocks_traj170_const_170_10_ini_mode_1.mat';
%
% test traj10_const_10_90_10_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_90_random_traj10_const_10_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_case2_traj10_const_10_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_SPC_ok_traj10_const_10_90_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_two_blocks_traj10_const_10_90_10_ini_mode_1.mat';
%
% test traj10_const_10_170_10_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_170_random_traj10_const_10_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_case2_traj10_const_10_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_SPC_ok_traj10_const_10_170_10_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_two_blocks_traj10_const_10_170_10_ini_mode_1.mat';
%
% test traj90_const_90_10_90_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_90_random_traj90_const_90_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_case2_traj90_const_90_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_SPC_ok_traj90_const_90_10_90_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_90_two_blocks_traj90_const_90_10_90_ini_mode_1.mat';
%
% test traj170_const_170_10_170_ini_mode_1
%
% string_data_sim ='Control_H_Data_square_10_170_random_traj170_const_170_10_170_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_case2_traj170_const_170_10_170_ini_mode_1.mat';
% string_data_sim ='Control_H_Data_square_10_170_SPC_ok_traj170_const_170_10_170_ini_mode_1.mat';
string_data_sim ='Control_H_Data_square_10_170_two_blocks_traj170_const_170_10_170_ini_mode_1.mat';

load([dataset_dir,string_data_sim])


%% Plots - Data Acquisition & Data structure Analysis

% Plotting the dataset
str_dataset_x = 't [s]';
str_dataset_u = '$u_d(t)$ [v]';
str_dataset_y = '$\theta_d(t)$ $[^\circ]$';
fig_data_traj= plotting_dataset(time_model, Data_u_model, 180/pi*Data_y_model, str_dataset_x, str_dataset_u, str_dataset_y, plt_set, 1);



% Plotting the dataset highlighting data belonging regions
plotting_dataset_hl(time_model, modes, Data_u_model, 180/pi*Data_y_model, str_dataset_x, str_dataset_u, str_dataset_y, plt_set, 2);

next_fig_index = 3;
fig_map = Hankel_map(H_y, modes, next_fig_index, 'Hankel');

% Identify the submatrixes for each mode
% sub_ind_H: matrix containing the indexes in the hankel matrix of the submatrices (if any)
% sub_num: number of submatrices
% fig_H_sub: figure handle of the submatrices plot
[~, sub_ind_H, sub_num, fig_H_sub] = Hankel_sub(H_y, modes, N_h, n_x, n_u, fig_map, 4);

% Adding the points to track to the map of the data
% count_data_N_h: counter of the number of data in N_h belonging to each mode or region amongst them (directionality in the switching is accounted) 
Hankel_Map_Append_id(H_y, [y_r_ini y_r_f(1:L)'], modes, N_ini, fig_map);


%% Plots - Reference

%Reference Input
fig_input = figure(5);
plot(time, u_r_f(1:L),'--k', 'linewidth',plt_set.thick), hold on, grid on,
% legend('Actual', 'Location','best'); grid on
x_title = xlabel(str_xaxis,'Interpreter','latex');
str_y = '$u_r(t)$ [v]';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Reference input','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_input.Units               = 'centimeters';
fig_input.Position(3)         = 9.5;
fig_input.Position(4)         = 7;
set(fig_input.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_input.PaperPositionMode   = 'auto';

% Reference Output
fig_output = figure(6);
plot(time, 180/pi*y_r_f(1:L),'--k','linewidth', plt_set.thick), hold on, grid on,
x_title = xlabel(str_xaxis,'Interpreter','latex');
str_y = '$\theta_r(t)$ $[^\circ]$';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Reference Output','Interpreter','latex', 'FontWeight', 'bold'); 
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_output.Units               = 'centimeters';
fig_output.Position(3)         = 9.5;
fig_output.Position(4)         = 7;
set(fig_output.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_output.PaperPositionMode   = 'auto';

%% Plots - Performances

% RMSE_u input plot (control period)
title_RMSE_u = 'System input $RMSE_u$ [v] (control period)';
x_axis_RMSE_u = '$\lambda_g$';
y_axis_RMSE_u = ['$RMSE_u$(',x_axis_RMSE_u,')'];
[RMSE_u_min, RMSE_u_max, lambda_g_RMSE_u_min, lambda_g_RMSE_u_max, ~] = Metric_Graph(lambda_g, RMSE_u, title_RMSE_u, x_axis_RMSE_u, y_axis_RMSE_u, 7);
% RMSE_y output plot (control period)
title_RMSE_y = 'System output $RMSE_y$ [rad] (control period)';
x_axis_RMSE_y = '$\lambda_g$';
y_axis_RMSE_y = ['$RMSE_y$(',x_axis_RMSE_y,')'];
[RMSE_y_min, RMSE_y_max, lambda_g_RMSE_y_min, lambda_g_RMSE_y_max, ~] = Metric_Graph(lambda_g, RMSE_y, title_RMSE_y, x_axis_RMSE_y, y_axis_RMSE_y, 8);
% RMSE_s plot (first instant of prediction vs real output)
title_RMSE_s = 'Predicted output $RMSE$ [rad] (first predicted instant over the control period)';
x_axis_RMSE_s = '$\lambda_g$';
y_axis_RMSE_s = ['$RMSE_{s,\hat{y}}$( ',x_axis_RMSE_s,')'];
[RMSE_s_min, RMSE_s_max, lambda_g_RMSE_s_min, lambda_g_RMSE_s_max, ~] = Metric_Graph(lambda_g, RMSE_s, title_RMSE_s, x_axis_RMSE_s, y_axis_RMSE_s, 9);
% Number of not null columns (single instant)
% c_not_null_slc_inst: selected instant of the Control period for the not null columns analysis (where t=0 corresponds to 1)
c_not_null_slc_inst = 1;
title_c_not_null_inst = ['Number of not null columns in g at instant ', num2str(c_not_null_slc_inst)];
x_axis_c_not_null = '$\lambda_g$';
y_axis_c_not_null = 'Number';
[c_not_null_min_inst, c_not_null_max_inst, lambda_g_c_not_null_min_inst, lambda_g_c_not_null_max_inst, ~] = Metric_Graph(lambda_g, c_not_null(c_not_null_slc_inst,:), title_c_not_null_inst, x_axis_c_not_null, y_axis_c_not_null, 10);
% Number of not null columns (all Control period)
title_c_not_null_glob = 'Number of not null columns in g - Control period';
[c_not_null_min_glob, c_not_null_max_glob, lambda_g_c_not_null_min_glob, lambda_g_c_not_null_max_glob] = Metric_HeatMap_L(lambda_g, L, c_not_null, title_c_not_null_glob, 11);
% Number of not null columns (average over the Control period instants)
c_not_null_avg = mean(c_not_null,1);
title_c_not_null_avg_inst = 'Number of not null columns in g - Average over the Control period';
x_axis_c_not_null = '$\lambda_g$';
y_axis_c_not_null = 'Number';
[c_not_null_avg_min, c_not_null_avg_max, lambda_g_c_not_null_avg_min, lambda_g_c_not_null_avg_max, ~] = Metric_Graph(lambda_g, c_not_null_avg, title_c_not_null_avg_inst, x_axis_c_not_null, y_axis_c_not_null, 12);

if sub_num > 0

    % % Rank of the sparsified matrix (each subamtrix)
    % title_r_s = 'Rank of the sparsified matrix';
    % x_axis_r_s = '$\lambda_g$';
    % y_axis_r_s = 'Rank';
    % Metric_Graph_Rank_sub(lambda_g, r_s, r_b, sub_ind_H(:,3), sub_num, title_r_s, x_axis_r_s, y_axis_r_s, 7);
  
    % r_s_slc_inst: selected instant of the Control period for the rank of the sparsified matrix analysis (where t=0 corresponds to 1)
    r_s_slc_inst = 1;
    % Rank of the sparsified matrix (average in a certain instant)
    r_s_avg_inst = NaN*ones(n_h, 1);
    % Rank of the sparsified matrix (average over all instants)
    r_s_avg = NaN*ones(n_h, 1);

    for i = 1 : n_h
        r_s_avg_inst(i,1) = mean(r_s{i}(r_s_slc_inst,:));
        r_s_avg(i,1) =  mean(mean(r_s{i}));
    end
    % Rank of the sparsified matrix (average in a certain instant)
    title_r_s_avg_inst = ['Rank of the sparsfied matrix - Average at instant ', num2str(c_not_null_slc_inst)];
    x_axis_r_s_avg_inst = '$\lambda_g$';
    y_axis_r_s_avg_inst = 'Rank';
    [r_s_avg_inst_min, r_s_avg_inst_max, lambda_g_r_s_avg_inst_min, lambda_g_r_s_avg_inst_max, fig_r_s_avg_inst] = Metric_Graph(lambda_g, r_s_avg_inst, title_r_s_avg_inst, x_axis_r_s_avg_inst, y_axis_r_s_avg_inst, 13);
    % Rank of the sparsified matrix (average over all instants)
    title_r_s_avg = 'Rank of the sparsfied matrix - Average over the Control period';
    x_axis_r_s_avg = '$\lambda_g$';
    y_axis_r_s_avg = 'Rank';
    [r_s_avg_min, r_s_avg_max, lambda_g_r_s_avg_min, lambda_g_r_s_avg_max, fig_r_s_avg] = Metric_Graph(lambda_g, r_s_avg, title_r_s_avg, x_axis_r_s_avg, y_axis_r_s_avg, 14);
 
    figure(fig_r_s_avg);
    yline(r_b, "--", 'DisplayName', 'Rank for the behavioural id', 'Color', [0 0 0]);
    % Distance of the (average) rank to the rank for the behavioural id
    r_dist = abs(r_s_avg - r_b*ones(length(lambda_g),1));
    title_r_dist = 'Distance of the average rank from the behavioural identification one - Average over the Control period';
    x_axis_r_dist = '$\lambda_g$';
    y_axis_r_dist = 'Distance - Rank value';
    [r_dist_min, r_dist_max, lambda_g_r_dist_min, lambda_g_r_dist_max, ~] = Metric_Graph(lambda_g, r_dist, title_r_dist, x_axis_r_dist, y_axis_r_dist, 15);
else 
    r_dist = NaN;
end

% Plotting the metrics together
title_tgt = 'Hyperparameter - Hankel';
x_axis_tgt = '$\lambda_g$';
y_axis_tgt = '$RMSE_u$';
z_axis_tgt = '$RMSE_y$';
fig_RMSE_comb = dim3plotting(lambda_g, RMSE_u, RMSE_y, title_tgt, x_axis_tgt, y_axis_tgt, z_axis_tgt, 16);
str_RMSE_u_min = 'RMSE_u min';
str_RMSE_y_min = 'RMSE_y min';
minima_plt_combined(fig_RMSE_comb, lambda_g, RMSE_u, RMSE_y, lambda_g_RMSE_u_min(1,1), lambda_g_RMSE_y_min(1,1), str_RMSE_u_min, str_RMSE_y_min);

writetable(table(RMSE_u_min, lambda_g_ins_1_min, RMSE_y_min, lambda_g_ins_2_min, RMSE_s_min, lambda_g_ins_3_min),'Results.txt');
open('Results.txt');

%% Plots - Control Deployment

% Input Constraint
u_min = -0.25;
u_max = 0.25;

% Input - min tracking error
fig_input = figure(17);
plot(time, u_r_f(1:L),'--k', 'linewidth',plt_set.thick), hold on,
plot(time, u{index_lambda_g_min_1, 1}(N_ini+1:N_ini+L),'Color', "#54278f",'linewidth',plt_set.thick), grid on,
plot(time, u{index_lambda_g_min_2, 1}(N_ini+1:N_ini+L),'Color', "#756bb1",'linewidth',plt_set.thick),
plot(time, u{index_lambda_g_min_3, 1}(N_ini+1:N_ini+L),'Color', "#9e9ac8",'linewidth',plt_set.thick),
yline(u_min, 'Color', "#808080", 'HandleVisibility','off'),
yline(u_max, 'Color', "#808080", 'HandleVisibility','off'),
legend('Reference Input', ['$RMSE_u$ min - $\lambda_g$ =', num2str(lambda_g_ins_1_min)],['$RMSE_y$ min - $\lambda_g$ =', num2str(lambda_g_ins_2_min)], ['$RMSE_{s,\hat{y}}$ min - $\lambda_g$ =', num2str(lambda_g_ins_3_min)], 'Location','best','Interpreter','latex'); grid on
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
str_y = '$u(t)$ [v]';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Input tracking','Interpreter','latex', 'FontWeight', 'bold'); % Salva l'handle del titolo
set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_input.Units               = plt_set.plot_unit;
fig_input.Position(3)         = plt_set.plot_dim_x;
fig_input.Position(4)         = plt_set.plot_dim_y;
set(fig_input.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_input.PaperPositionMode   = 'auto';

% Output - min tracking error
fig_output = figure(18);
plot(time, 180/pi*y_r_f(1:L),'--k','linewidth',plt_set.thick), hold on,
plot(time, 180/pi*y{index_lambda_g_min_1, 1}(N_ini+1 : end-1, 1),'Color', "#54278f",'linewidth',plt_set.thick), 
plot(time, 180/pi*y{index_lambda_g_min_2, 1}(N_ini+1 : end-1, 1),'Color', "#756bb1",'linewidth',plt_set.thick), 
plot(time, 180/pi*y{index_lambda_g_min_3, 1}(N_ini+1 : end-1, 1),'Color', "#9e9ac8",'linewidth',plt_set.thick), 
legend('Reference Output', ['$RMSE_u$ min - $\lambda_g$ =', num2str(lambda_g_ins_1_min)],['$RMSE_y$ min - $\lambda_g$ =', num2str(lambda_g_ins_2_min)], ['$RMSE_{s,\hat{y}}$ min - $\lambda_g$ =', num2str(lambda_g_ins_3_min)], 'Location','best','Interpreter','latex'); grid on
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
str_y = '$\theta(t)$ $[^\circ]$';
y_title = ylabel(str_y,'Interpreter','latex');
str_title = 'Output tracking';
h_title = title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');set(h_title, 'FontSize', plt_set.title_font_dim);
set(y_title, 'FontSize', plt_set.y_font_dim);
set(x_title, 'FontSize', plt_set.x_font_dim);
fig_output.Units               = plt_set.plot_unit;
fig_output.Position(3)         = plt_set.plot_dim_x;
fig_output.Position(4)         = plt_set.plot_dim_y;
set(fig_output.Children, ...
    'FontName',     plt_set.fontname, ...
    'FontSize',     plt_set.fontsize);
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02))
fig_output.PaperPositionMode   = 'auto';

% Plotting the weights of the trajectories
[fig_weights_all, fig_weights_avg_var] = plotting_weights_L(g_opt{index_lambda_g_min_2, 1}, lambda_g(index_lambda_g_min_2), plt_set, 19, 20);

%% Plots - Predictions

% index_lambda_g_pre: index of lambda_g to inspect
% (index_lambda_g_min _1: min RMSE_u, _2: min RMSE_y, _3: min RMSE_s)
index_lambda_g_pre = index_lambda_g_min_2;
% Color_re_pre: plot color of the real output plotted alongside the
% predictions (It should correspond with index_lambda_g_pre)
Color_re_pre = "#756bb1";
% string_re_pre: string of the real output plotted alongside the
% predictions (It should correspond with index_lambda_g_pre)
string_re_pre = ['$RMSE_y$ min - $\lambda_g$ =', num2str(lambda_g_ins_2_min)];

% inst: Instant to inspect of the Control period
% We choose the initial instant, the one just before the switch and the one just after it HERE
% one switch
% inst = [0 14 15]';
% % two switches CONTROLLARE CHE SONO GIUSTI
inst = [0 9 10 19 20]'; 


str_u_pre_inst_y = '$u(t)$';
str_u_pre_inst_var_name = 'Input';
fig_u_pre_inst = predict_inst_plt_var(inst, time, u_opt_seq{index_lambda_g_pre, 1}, u{index_lambda_g_pre, 1}(N_ini+1:N_ini+L), u_r_f, L, tau_s, lambda_g(index_lambda_g_pre), str_xaxis, str_u_pre_inst_y, str_u_pre_inst_var_name, plt_set, Color_re_pre, string_re_pre, 21);
% Adding constraints
figure(fig_u_pre_inst)
yline(u_min, 'Color', "#808080", 'HandleVisibility','off'),
yline(u_max, 'Color', "#808080", 'HandleVisibility','off'),
str_y_pre_inst = '$y(t)$';
str_y_pre_inst_var_name = 'Output';
fig_y_pre_inst = predict_inst_plt_var(inst, time, 180/pi*y_pre_seq{index_lambda_g_pre, 1}, 180/pi*y{index_lambda_g_pre, 1}(N_ini+1 : end-1, 1), 180/pi*y_r_f, L, tau_s, lambda_g(index_lambda_g_pre), str_xaxis, str_y_pre_inst, str_y_pre_inst_var_name, plt_set, Color_re_pre, string_re_pre, 22);
str_dec_x = 'Column';
str_dec_y = 'Weight';
str_dec_var_name = 'Weights';
fig_g_opt_inst = predict_inst_plt_dec(inst, g_opt{index_lambda_g_pre, 1}, lambda_g(index_lambda_g_pre), str_dec_x, str_dec_y, str_dec_var_name, 23);