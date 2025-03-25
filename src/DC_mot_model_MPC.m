function zdot=DC_mot_model_MPC(z,u,e, tau_s, th)

%% Read parameters, states and inputs

% Parameters
g       =       th(1,1);     %
J       =       th(2,1);     % 
m       =       th(3,1);     % 
k_m     =       th(4,1);     % 
l       =       th(5,1);     % 
tau     =      th(6,1);     % 

% States
theta_start               =       z(1,1);    % [rad]
theta_dot_start           =       z(2,1);    % [rad/s]

% Model equations

theta_dot = theta_dot_start;
theta_dotdot =  (m*g*l)/J * cos(theta_start) - 1/tau*theta_dot_start + k_m/tau * u;

% theta_dotdot =  (m*g*l)/J * (4/pi^2*theta_start^2-8/pi*theta_start+3) - 1/tau*theta_dot_start + k_m/tau * u;

% Discretization via Forward Euler
    zdot = [ theta_start + tau_s * theta_dot;
             theta_dot_start + tau_s * theta_dotdot];







