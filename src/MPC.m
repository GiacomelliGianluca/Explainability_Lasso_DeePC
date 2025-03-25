function [u_seq_opt, y_seq_opt, J_opt] = MPC(u_r, y_r, x_init, th, n_u, n_y, n_x, Q, R, N_h, tau_s, N_ex)
  
% FHCOP_rate
    % Solve the Finite Horizon Control Optimization Problem
    
    % Parameters
    % - u_r: The target input
    % - y_r: The target output
    % - x_init: initial state
    % - th: system parameters
    % - n_u: number of input
    % - n_y: number of output
    % - n_x: number of states
    % - Q: The state weight
    % - R: The input weight
    % - N_h: The prediction horizon (number of steps)
    % - tau_s: The sampling time
    % - N_ex: Execution step of MPC
    % - time: current time
    % - N_i: current step (overall)

    % x_fin:terminal constraints
    x_fin = [y_r(end);0];

    % Import the CasADi toolkit
    import casadi.*;
    opti = casadi.Opti();
    
   
    % Declare the optimization variables
    u = opti.variable(n_u,  N_h);
    x = opti.variable(n_x, N_h + 1);

    % Cost function initialization
    J = 0;

    % Initial state
    opti.subject_to(x(:, 1) == x_init);
    % Input Constraints
    % opti.subject_to(-0.25 <= u <= 0.25);
    % opti.subject_to(-5 <= u <= -0.5);
    % opti.subject_to(0.5 <= u <= 5);
    opti.subject_to(-10 <= u <= 10);


    opti.subject_to(x(:, N_h + 1) == x_fin);

            
    for ii=1:N_h

        % Compute x(k+1) = x(k) + f(x(k), u(k)) * Ts
        xp = DC_mot_model_MPC(x(:,ii),u(1,ii),0, tau_s, th);

        opti.subject_to((x(:, ii+1) - xp) == 0);


                % Cost function 
            J = J ...
                + (u(1,((ii-1)*n_u+1):ii*n_u) - u_r(ii)).' * R * (u(1, (ii-1)*n_u+1:ii*n_u) - u_r(ii)) ...
                + (x(1, ii) - y_r(ii)).' * Q * (x(1, ii) - y_r(ii));

          % opti.subject_to(-5 <= u(ii) - u(ii+1) <= 5);
    end



    % Set the initial guess
    % x(k+i|k) = x(k) ∀i = 1, ..., N
    opti.set_initial(x, repmat(x_init, 1, N_h+1));    
    opti.set_initial(u, u_r);   

    % Declare the cost function
    opti.minimize(J);                   

    % CASADI settings: da fare se vuoi togliere il printing
          prob_opts = struct;
    %     prob_opts.expand = true;
          prob_opts.ipopt.print_level = 0;    % Disable printing 0
          prob_opts.print_time = false;       % Do not print the timestamp
            
    
    %%%%%% CasADi Settings (do not change) %%%%%%
   
    % Options of the solver
    ip_opts = struct;
    ip_opts.print_level = 0;            % Disable printing
    ip_opts.max_iter = 1e5;             % Maximum iterations
    ip_opts.compl_inf_tol = 1e-6;
    
    % Set the solver
    opti.solver('ipopt', prob_opts, ip_opts);

    try
        % SOLVE THE FHOCP
        sol = opti.solve();

        % Extract the optimal control action
        u_seq_opt = sol.value(u(1, 1:n_u*N_ex));
        y_seq_opt = sol.value(x(1, 1:n_y*N_ex));
        J_opt = sol.value(J);


    catch EX
        opti.debug.show_infeasibilities() %In case the solver stops due to problem infeasibility, you may identify the problematic constraints with:
        keyboard;
    end



end