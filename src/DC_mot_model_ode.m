function dydt=DC_mot_model_ode(t,y,u,e,th)

%% Read parameters, states and inputs

% Parameters
g       =       th(1,1);     %
J       =       th(2,1);     % 
m       =       th(3,1);     % 
k_m     =       th(4,1);     % 
l       =       th(5,1);     % 
tau     =      th(6,1);     % 

% States
% y(1): theta [rad]
% y(2): theta_dot [rad/s]

% Model equations
dydt = [y(2);
        (m*g*l)/J * cos(y(1)) - 1/tau*y(2) + k_m/tau * u ];





