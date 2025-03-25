% Main script for the unbalnaced disk benchmark of 
% "Insights into the explainability of Lasso-based DeePC for nonlinear
% systems", G. Giacomelli, S. Formentin, V. Lopez, M. Müller, and V. Breschi

close all 
clear 
clc

% plt_set: Plot settings
plt_set.y_font_dim = 16;
plt_set.x_font_dim = 11;
plt_set.title_font_dim = 14;
plt_set.thick = 1.6;
plt_set.plot_dim_x = 8; 
plt_set.plot_dim_y = 6; 
plt_set.plot_unit = 'centimeters'; 
plt_set.fontname = 'Times';
plt_set.fontsize = 12;

% Dataset directory
dataset_dir = 'Please insert data directory';

%% System %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n_x number of states
n_x = 2;
% n_u Input Size
n_u = 1;
% n_y Output Size
n_y = 1;
 

% Parameters of the unbalanced disk 
g       =       9.81;     % Gravitational acceleration [m/s^2]
J       =       2.4e-04;  % Disk inertia [Nm^2]  
m       =       0.076;    % Mass [kg]
k_m     =       10.5;     % Motor constant [-]    
l       =       0.041;    % Distance [m]
tau     =       0.4;      % Lumped back EMF constant
% Parameter vector
th       =       [g;J;m;k_m;l;tau];

%% Reading the data

% Selecting the datasest to build the Hankel 

% string_data_hankel = 'Data_017.mat';
% string_data_hankel = 'Data_017_157_multiple_transitions.mat';
% string_data_hankel = 'Data_017_157_single_transition.mat';
string_data_hankel = 'Data_017_297_multiple_transitions.mat';
% string_data_hankel = 'Data_017_297_single_transition.mat';
Data_u_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "Data_u")));
Data_y_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "Data_y")));
time_model = cell2mat(struct2cell(load([dataset_dir,string_data_hankel], "time")));

% Dataset length (Assuming that Data_u and Data_y have the same data)
N_D = length(Data_u_model);

%% Partition the data

% Forming the data vector
Data = [Data_u_model' Data_y_model'];

% Plotting the data distribution
% Declaring the Operating Points (in deg) 
OPs = [10;
       170];

% Declaring the condition
cond.y = mean(OPs);

% Initial data plot
Intial_data = figure(1);
plot(Data(:,1),Data(:,2),'k*','MarkerSize',3);
grid on,
str_x = '$u_d(t)$ [V]';
xlabel(str_x,'Interpreter','latex');
str_y = '$y_d(t)$ [rad]';
ylabel(str_y,'Interpreter','latex');
str_title = 'Dataset';
title(str_title, 'Interpreter','latex', 'FontWeight', 'bold');

% n_p: Number of operating points
% (set to 1 for the experiment with a single OP and comment coherently the ensuing part)
n_p = 2;

% Partitioning data
[label, Centr] = Partitioning(Data_y_model, cond.y*pi/180);


% Partitioned data plot 
title_string = 'Partitioned dataset';
Clustered_data_id = Partitioned_data_plot(Data, label, title_string, str_x, str_y, 2);

%% DDPC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters of DDPC %

% Past horizon
N_ini = 40;
% Future horizon 
L  = 30;
% Prediction Horizon
N_h = N_ini + L;

% Control period
N_c = 100;

% Sampling time [s] 
tau_s = 1e-2; 


% Number of simulation instants
N = N_ini + N_c; 

% Execution period of DDPC [s]
t_ex = tau_s;
% Number of steps where the DDPC is not updated
N_ex = t_ex/tau_s;

% Setting the references [rad] 
% Output reference evolving from x_0 to x_final based on spline % Initial point of the trajectory
y_0_r_ini_deg = 10;
y_0_r_ini = pi/180*y_0_r_ini_deg;
% Final point of the trajectory
y_f_r_ini_deg = 10;
y_f_r_ini = pi/180*y_f_r_ini_deg;
% Defining the middle middle points for spline (spline since it is cubic
% needs 4 points to work)
if y_0_r_ini < y_f_r_ini
    y_mid_1_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.25 + y_0_r_ini;
    y_mid_2_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.75 + y_0_r_ini;
else
    y_mid_1_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.75 + y_f_r_ini;
    y_mid_2_r_ini = abs(y_0_r_ini - y_f_r_ini)*0.25 + y_f_r_ini;
end
% Defining the points for spline
y_spline = [ y_0_r_ini  y_mid_1_r_ini  y_mid_2_r_ini   y_f_r_ini];
x_spline = [1 N_ini/4  3/4*N_ini N_ini];
% Defining the reference vector
y_r_ini = spline(x_spline, y_spline, 1:1:N_ini);

% % Input reference (past horizon steps) (unused, see R_MPC)
u_r_ini = computing_u_equil_pen(y_r_ini, th);

% Periods of the reference signal 
N_r_1 = round(L/2); % Before the switch 
N_r_2 = N_c - N_r_1; % After the switch 

% Future horizon references 
% Output reference sequence  
y_r_1_deg = 10*ones(1,N_r_1);
y_r_1 = pi/180*y_r_1_deg;
y_r_2_deg = 170*ones(1, N_r_2);
y_r_2 = pi/180*y_r_2_deg;
% y_r_ch: output reference control horizon
y_r_ch = [y_r_1     y_r_2];
% Checking the correct length
if length(y_r_ch)~=N_c
    keyboard;
end
y_r_f = [y_r_ch  y_r_ch(end)*ones(1,L)]';

% Input reference
u_r_f = computing_u_equil_pen(y_r_f, th);

% Time vector [s] (x-axis)
time = 0: tau_s: (N_c - 1)*tau_s;

str_xaxis = '$t$ [s]';

%% MPC's weights for the initial trajectory
% Q Matrix (weight of y) for MPC tracking of the past horizon traj
Q_MPC = 1e3*eye(n_y);
% R Matrix (weight of u) for MPC tracking of the past horizon traj
R_MPC = 0*eye(n_u);

%% DeePC's weights and hyperparameter
% Q Matrix (weight of y) for DeePC
Q_DeePC = 1e2*eye(n_y);
% R Matrix (weight of u) for DeePC
R_DeePC = 1*eye(n_u);

%Variation for the hyperparameter
n_h = 11; 
lambda_g = logspace(-5, 5, n_h);

% Input Hankel Matrix
H_u = building_Hankel(Data_u_model, N_h, n_u, n_x, 1);
% Output Hankel Marix
H_y = building_Hankel(Data_y_model, N_h, n_y, n_x, 0);
% Hankel with the labels (SISO case)
H_label = building_Hankel(label, N_h, 1, 1, 0);

next_fig_index = 3;
fig_map = Hankel_map(OPs, H_label, plt_set, next_fig_index, 'Hankel');

% Identify the submatrixes for each OP
% sub_ind_H: matrix containing the indexes in the hankel matrix of the submatrices (if any)
% sub_num: number of submatrices
% fig_H_sub: figure handle of the submatrices plot
[~, sub_ind_H, sub_num, fig_H_sub] = Hankel_sub(H_y, H_u, OPs, N_h, n_x, n_u, cond, fig_map, plt_set, 4);

% Allocation input applied (First N_ini steps is determined by an MPC, then
% the remaining L steps is decided through DDPC)
u = cell(n_h, 1);
% Allocation system output (First N_ini steps is determined by an MPC, then
% the remaining L steps is decided through DDPC)
y = cell(n_h, 1);
% Allocation overall cost
J_opt = cell(n_h, 1);
% Allocation tracking cost
J_tr = cell(n_h, 1);
% Allocation reg on g cost
J_g = cell(n_h, 1);

% Allocation of the optimsl inpute sequence (computed by the DDPC)
u_opt_seq = cell(n_h, 1);
% Allocation of the predicted output (computed by the DDPC)
y_pre_seq = cell(n_h, 1);
% Allocation of optimal weights 
g_opt = cell(n_h, 1);

% Initializing the predicted output and weights
for k = 1 :n_h
    % For every element of the cell is allocated the matrix corresponding to the optimal input applied or the consequent system output for the different lambda_g_k 
    u{k,1} = NaN*ones(n_u*(N_ini+N_c),1);
    y{k,1} = NaN*ones(n_y*(N_ini+N_c),1);
    % Similarly for the costs and the metrics on g
    J_opt{k,1} = NaN*ones(N_c,1);
    J_tr{k,1} = NaN*ones(N_c,1);
    J_g{k,1} = NaN*ones(N_c,1);
    % For every element of the cell is allocated the matrix corresponding to the optimal input sequence or predicted output for the different lambda_g_k 
    u_opt_seq{k,1} = NaN*ones(n_u*L,N_c);
    y_pre_seq{k,1} = NaN*ones(n_y*L,N_c);
    g_opt{k,1} = NaN*ones(N_D-N_h+1,N_c);
end


%% Metrics  
% RMSE_u (Prediction horizon) 
RMSE_u = NaN*ones(length(lambda_g), 1);
% RMSE_y (Prediction horizon) 
RMSE_y = NaN*ones(length(lambda_g), 1);

%% Execution

% Allocation of the vector of the simulation (each column corresponds to a state)
y_sim = zeros(N, n_x);
% Initialization IC ode_45
y_sim (1, 1) = y_r_ini(1,1);

for i=1:length(lambda_g)
    % Reallocation of u_ini, y_ini
    u_ini = zeros(n_u*N_ini, 1);
    y_ini = NaN*ones(n_y*N_ini, 1);
    for k= 1 : (N_ini + N_c)
    
        % Past horizon
        if k <= N_ini
            % MPC case 
            if k == 1
                % MPC for the trajectory of the past horizon (prediction horizon N_ini)
                [u_c_seq, ~, ~] = MPC(u_r_ini, y_r_ini, y_sim(k, :)', th, n_u, n_y, n_x, Q_MPC, R_MPC, N_ini, tau_s, N_ini);
             end
             % Update MPC control actions
             u{i,1}(k) = u_c_seq(1,((k-1)*n_u+1):n_u*k);
         end
    
        
         % Computation of the optimal control action
         if k >= N_ini + 1          % Checking if the DeePC has to be executed
             % Index control period
             fut = k - N_ini;
             % Solve the DeePC
             % Lasso-DeePC
             [J_opt{i,1}(fut), J_tr{i,1}(fut), J_g{i,1}(fut), g_opt{i,1}(:,fut), u_opt_seq{i,1}(:,fut), y_pre_seq{i,1}(:,fut)] = DeePC_y_f_lambda_g(u_ini, y_ini, u_r_f(fut: fut + L - 1), y_r_f(fut: fut + L - 1), n_u, n_y, Q_DeePC, R_DeePC, N_h, N_ini, N_ex, H_u, H_y, N_D,  lambda_g(i));
             % Update control action
             u{i,1}(k) = u_opt_seq{i,1}(1,fut);
         end  
         
         % System evolution
         [~,y_ode] = ode45(@(t,y) DC_mot_model_ode(t,y, u{i,1}(k),0, th), [0 tau_s], y_sim (k, :));
        
         % Update of the initial traj
         u_ini = [u_ini(n_u*2:end,1);    u{i,1}(k)];
         y_ini = [y_ini(n_u*2:end,1);    y_sim(k, 1)];
    
         y_sim (k+1, :) = y_ode(end, :);
            
     end
    % Saving the simulation output
    y{i,1} = y_sim;
    
    if~isnan(y_pre_seq{i,1})
       % Metrics computation 
       % RMSE_u input (Prediction horizon) 
       RMSE_u(i) = sqrt(1/N_c*sum((u{i,1}(N_ini+1 : end, 1) - u_r_f(1: N_c)).^2));        
       % RMSE_y output (Prediction horizon) 
       RMSE_y(i) = sqrt(1/N_c*sum((y{i,1}(N_ini+1 : end-1, 1) - y_r_f(1: N_c)).^2));
    end

end

%% Metric Plotting

% RMSE_u input plot (control period)
title_RMSE_u = 'System input $RMSE_u$ [V] (control period)';
x_axis_RMSE_u = '$\lambda_g$';
y_axis_RMSE_u = ['$RMSE_u$(',x_axis_RMSE_u,')'];
[RMSE_u_min, RMSE_u_max, lambda_g_RMSE_u_min, lambda_g_RMSE_u_max, ~] = Metric_Graph(lambda_g, RMSE_u, title_RMSE_u, x_axis_RMSE_u, y_axis_RMSE_u, 0);
% RMSE_y output plot (control period)
title_RMSE_y = 'System output $RMSE_y$ [rad] (control period)';
x_axis_RMSE_y = '$\lambda_g$';
y_axis_RMSE_y = ['$RMSE_y$(',x_axis_RMSE_y,')'];
[RMSE_y_min, RMSE_y_max, lambda_g_RMSE_y_min, lambda_g_RMSE_y_max, ~] = Metric_Graph(lambda_g, RMSE_y, title_RMSE_y, x_axis_RMSE_y, y_axis_RMSE_y, 0);

title_tgt = 'Hyperparameter - Hankel';
x_axis_tgt = '$\lambda_g$';
y_axis_tgt = '$RMSE_u$';
z_axis_tgt = '$RMSE_y$';
fig_RMSE_comb = dim3plotting(lambda_g, RMSE_u, RMSE_y, title_tgt, x_axis_tgt, y_axis_tgt, z_axis_tgt, plt_set, 5);

%% Plotting the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The results will be plotted for the selected lambda_g (In case of
% multiple, the first one will be considered)
% Minimum RMSE_y 
lambda_g_ins_min = lambda_g_RMSE_y_min(1,1);
% Retrieving the indexes
index_lambda_g_min = find(lambda_g == lambda_g_ins_min);

% Input Constraint
u_min = -10;
u_max = 10;

% Input
fig_input = figure(6);
plot(time, u_r_f(1:N_c),'-k', 'linewidth',1), hold on, grid on,
plot(time, u{index_lambda_g_min, 1}(N_ini+1:N_ini+N_c),'Color', "#9e9ac8",'linewidth',plt_set.thick),
yline(u_min, 'Color', [0.5 0.5 0.5], 'linewidth',1.5, 'LineStyle',':' ,'HandleVisibility','off'),
yline(u_max, 'Color', [0.5 0.5 0.5], 'linewidth',1.5, 'LineStyle',':' ,'HandleVisibility','off'),
ylim([-11.75 11.75])
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
str_y = '$u(t)$ [V]';
y_title = ylabel(str_y,'Interpreter','latex');
h_title = title('Input Tracking','Interpreter','latex', 'FontWeight', 'bold');
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

% Output 
fig_output = figure(7);
plot(time, y_r_f(1:N_c),'-k','linewidth',1), hold on, grid on,
plot(time, y{index_lambda_g_min, 1}(N_ini+1 : end-1, 1),'Color', "#9e9ac8",'linewidth',plt_set.thick), 
str = str_xaxis;
x_title = xlabel(str,'Interpreter','latex');
str_y = '$y(t)$ [rad]';
y_title = ylabel(str_y,'Interpreter','latex');
str_title = 'Output Tracking';
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
% % inst: Instant to inspect of the control period
% inst = [0 15]';
inst = [0 60]';
% 
str_dec_x = '$k$ [column]';
str_dec_y = '$[g]_k$';
str_dec_var_name = 'Data Selection';
fig_g_inst = predict_inst_plt_dec_p(inst, g_opt{index_lambda_g_min, 1}, lambda_g(index_lambda_g_min), sub_num, sub_ind_H, str_dec_x, str_dec_y, str_dec_var_name, plt_set, 8);

%% Saving performances in txt
writetable(table(RMSE_u_min,RMSE_y_min, lambda_g_ins_min),'Results.txt');
open('Results.txt');