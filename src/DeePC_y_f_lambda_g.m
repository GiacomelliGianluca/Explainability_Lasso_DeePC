function [J_all, J_tr, J_g, g_opt, u_opt_seq, y_opt_seq] = DeePC_y_f_lambda_g(u_ini, y_ini, u_r, y_r, n_u, n_y, Q, R, N_h, N_ini, N_ex, H_u, H_y, D, lambda_g)

  
    % Solve the DeePC optimization problem using the mosaic Hankel matrixes
    
    % Parameters
    % - u_ini: previous input
    % - y_ini: previous output
    % - u_r: The target input 
    % - y_r: The target output
    % - Q: The state weight
    % - R: The input weight
    % - N_h: The prediction horizon (number of steps)
    % - N_ini: Number of initialization steps
    % - H_u: Hankel matrix input
    % - H_y: Hankel matrix output
    % - D: number of data


    % Future Horizon
    L_h = N_h-N_ini;


    % Start cvx
    cvx_begin quiet
        % Declare the optimization variables
        variable u(n_u*L_h)
        variable y(n_y*L_h)
        variable g(D-N_h+1)
        % variable sigma_y(n_y*N_ini,1)           

        % Cost function initialization
        J = 0;
        
        for ii=1:L_h

            % Cost function 
            J = J ...
                + (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)).' * Q * (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)) ...
                + (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1)).' * R * (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1));
        end
        % Adding the regulizer at the cost function
        J = J + lambda_g*norm(g,1);

        minimize(J)
        subject to
            % Model constraints
            H_u(1:N_ini,:)*g==u_ini;
            H_u(N_ini+1:end,:)*g==u;
            % H_y(1:N_ini,:)*g == y_ini + sigma_y;
            H_y(1:N_ini,:)*g == y_ini;
            H_y(N_ini+1:end,:)*g==y;
            % % Terminal constraints
            % opti.subject_to(u(end-1:end) == u_r);
            % opti.subject_to(y(end-1:end) == y_r);
            
            % Input constraints 
            -10 <= u <= 10


    cvx_end

    % Extract the cost function and the optimal weightings
    J_all = J;
    J_tr = 0;
    for ii=1:L_h
        % Cost function 
        J_tr = J_tr ...
             + (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)).' * Q * (y(((ii-1)*n_y+1):ii*n_y, 1) - y_r(((ii-1)*n_y+1):ii*n_y, 1)) ...
             + (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1)).' * R * (u(((ii-1)*n_u+1):ii*n_u, 1) - u_r(((ii-1)*n_u+1):ii*n_u, 1));
    end
    J_g = norm(g,1);
    % J_sig = norm(sigma_y,1);
    g_opt = g;
    % Extract the optimal input
    u_opt_seq = u;
    % Extract the optimal output prediction (for the sequence)
    y_opt_seq = y;

    if strcmp(cvx_status, 'Failed')
        disp('The problem is infeasible.')
        keyboard;
    end



end