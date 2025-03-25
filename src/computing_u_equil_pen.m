function u_bar = computing_u_equil_pen(y_bar, th)
% The function computes the corresponding input of a desired output of a
% PWA system

% y_bar: equilibrium outputs
% th: system parametes

% u_bar: equilibrium inputs

% Parameters
g       =       th(1,1);     %
J       =       th(2,1);     % 
m       =       th(3,1);     % 
k_m     =       th(4,1);     % 
l       =       th(5,1);     % 
tau     =      th(6,1);      %  

% Allocating u_bar
u_bar = NaN*ones(length(y_bar),1);

% Computing the equilibria
for i = 1 : length(y_bar)
    % Model equations
    u_bar(i) = -(tau/k_m)* (m*g*l/J) * cos(y_bar(i));
end 




end

