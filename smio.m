function [z_k_underline_trajectory, z_k_bar_trajectory] = ...
    smio(f, f_d, g, y, u, x_0_underline, ...
         x_0_bar, d_0_underline, d_0_bar, w_underline, w_bar, ...
         v_underline, v_bar, L_f, L_g, L_h, ...
         num_grid_points_per_dim_f_domain, ...
         num_grid_points_per_dim_g_domain, ...
         num_grid_points_per_dim_h_domain, K, ...
         num_interval_optimizations, varargin)
    % SMIO from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    %
    % Args:
    %   f: Function handle of a vector field f: R^{n + p + m + n_w} -> R^n.
    %   f_d: Function handle of the decomposition function of f,
    %        f_d: R^{n + p + m + n_w} x R^{n + p + m + n_w} -> R^n.
    %   g: Function handle of a vector field g: R^{n + p + m + l} -> R^l.
    %   y: A matrix of size (K + 1 x l) of obervations.
    %      The first entry is the observation of step 0.
    %   u: A matrix of size (K + 1 x m) of input vectors.
    %      The first entry is the input of step 0.
    %   L_f: Vector of Lipschitz constants for each codomain dimension of 
    %        f with shape (n x 1).
    %   L_g: Vector of Lipschitz constants for each codomain dimension of 
    %        g with shape (l x 1).
    %   L_h: Vector of Lipschitz constants for each codomain dimension of 
    %        h with shape (p x 1).
    %   K: The number of observations and corresponding iterations k.
    %   num_interval_optimizations: The number of optimizations i of the 
    %                               interval within each observation 
    %                               step k.
    % varargin:
    %   Key-Word Args:
    %       plot_global_affine_abstraction_f: Default false, type boolean.
    %       plot_global_affine_abstraction_g: Default false, type boolean.
    %       plot_global_affine_abstraction_h: Default false, type boolean.
    %       h: Function handle of a vector field 
    %          h: R^{n + p + m + n_w} -> R^p. Not necessary for smio. If
    %          provided, it will be plotted.
    %       x_spacing: Default 0.01. The spacing between plotting grid 
    %                  points on the x-axis.
    %       y_spacing: Default 0.01. The spacing between plotting grid 
    %                  points on the y-axis.
    %
    % Returns:
    %   - z_k_underline_trajectory: The calculate trajectory of lower 
    %                               bounds of z. Array of size 
    %                               (K + 1 x n + p).
    %   - z_k_bar_trajectory: The calculate trajectory of upper bounds of
    %                         z. Array of size (K + 1 x n + p).
    % Note: The bb in variable names is an abbreviation for blackboard 
    %       bold.
    %% Parse Inputs
    input_parser = inputParser;
    addParameter(input_parser, 'plot_global_affine_abstraction_f', ...
                 false, @islogical);
	addParameter(input_parser, 'plot_global_affine_abstraction_h', ...
                 false, @islogical);
	addParameter(input_parser, 'plot_global_affine_abstraction_g', ...
                 false, @islogical);
    addParameter(input_parser, 'h', -1);
    addParameter(input_parser, 'x_spacing', 0.01);
    addParameter(input_parser, 'y_spacing', 0.01);
    addParameter(input_parser, 'verbose', false, @islogical);
    addParameter(input_parser, 'log', @disp);
    parse(input_parser, varargin{:});
    results = input_parser.Results;
    plot_global_affine_abstraction_f = ...
        results.plot_global_affine_abstraction_f;
    plot_global_affine_abstraction_h = ...
        results.plot_global_affine_abstraction_h;
    plot_global_affine_abstraction_g = ...
        results.plot_global_affine_abstraction_g;
    h = results.h;
    x_spacing = results.x_spacing;
    y_spacing = results.y_spacing;
    verbose = results.verbose;
    log = results.log;
    
    %% Initialization
    n = size(x_0_underline, 1);
    m = size(u, 2);
    n_w = size(w_underline, 1);
    l = size(v_underline, 1);
    p = size(d_0_underline, 1);
    z_0_underline = [x_0_underline;
                     d_0_underline];
    z_0_bar = [x_0_bar;
               d_0_bar];
    xi_0_underline = [z_0_underline;
                      u(1, :)';
                      w_underline];
    xi_0_bar = [z_0_bar;
                u(1, :)';
                w_bar];
    d_underline = zeros(K + 1, p);  % d_k_underline storage
    d_bar = zeros(K + 1, p);  % d_k_bar storage
    xi_tilde = zeros(K + 1, n + p + m + n_w);  % xi_tilde storage
    xi_tilde(1, :) = (1/2) * (xi_0_bar + xi_0_underline);
    epsilon = zeros(K + 1, p);  % epsilon storage
    epsilon(1, :) = 2 * L_h * norm(xi_0_bar - xi_0_underline);
    z_k_underline_trajectory = zeros(K + 1, n + p);
    z_k_bar_trajectory = zeros(K + 1, n + p);
    
    %% Global Parallel Affine Abstraction for f
    % TODO: Abstraction matricies do not fit the given ones in the paper
    %       and e goes overboard, why?
    sigma_upper_bound_f = compute_sigma_upper_bound(L_f, ...
                                                  xi_0_underline, ...
                                                  xi_0_bar);
    [bb_A_f, e_bar_f, e_underline_f] = ...
        global_parallel_affine_abstraction( ...
            f, f, xi_0_underline, xi_0_bar, ...
            num_grid_points_per_dim_f_domain, sigma_upper_bound_f, n);
    
    if plot_global_affine_abstraction_f
        func_plot_global_affine_abstraction_f(f, n, p, m, n_w, ...
                                              x_0_underline, ...
                                              x_0_bar, ...
                                              d_0_underline, ...
                                              d_0_bar, ...
                                              u, ...
                                              w_underline, ...
                                              w_bar, ... 
                                              x_spacing, ...
                                              y_spacing, ...
                                              bb_A_f, ...
                                              e_underline_f, ...
                                              e_bar_f);
    end

    %% Global Parallel Affine Abstraction for h
    h_0_bar = @(xi) d_0_bar;
    h_0_underline = @(xi) d_0_underline;
    sigma_upper_bound_h = compute_sigma_upper_bound(L_h, ...
                                                    xi_0_underline, ...
                                                    xi_0_bar);
    [bb_A_h, e_bar_h, e_underline_h] = ...
        global_parallel_affine_abstraction( ...
            h_0_underline, h_0_bar, xi_0_underline, xi_0_bar, ...
            num_grid_points_per_dim_h_domain, sigma_upper_bound_h, p);
        
    if plot_global_affine_abstraction_h
        func_plot_global_affine_abstraction_h(h, h_0_underline, ...
                                              h_0_bar, ...
                                              n, p, m, n_w, ...
                                              x_0_underline, ...
                                              x_0_bar, ...
                                              d_0_underline, ...
                                              d_0_bar, ...
                                              u, ...
                                              w_underline, ...
                                              w_bar, ... 
                                              x_spacing, ...
                                              y_spacing, ...
                                              bb_A_h, ...
                                              e_underline_h, ...
                                              e_bar_h);
    end
        
    %% Global Parallel Affine Abstraction for g
    xi_0_underline_g = [z_0_underline;
                        u(1, :)';
                        v_underline];
    xi_0_bar_g = [z_0_bar;
                  u(1, :)';
                  v_bar];
    sigma_upper_bound_g = compute_sigma_upper_bound(L_g, ...
                                                    xi_0_underline_g, ...
                                                    xi_0_bar_g);
    [bb_A_g, e_bar_g, e_underline_g] = ...
        global_parallel_affine_abstraction( ...
            g, g, xi_0_underline_g, xi_0_bar_g, ...
            num_grid_points_per_dim_g_domain, sigma_upper_bound_g, l);
        
    if plot_global_affine_abstraction_g    
        func_plot_global_affine_abstraction_g(g, n, p, m, l, ...
                                              x_0_underline, ...
                                              x_0_bar, ...
                                              d_0_underline, ...
                                              d_0_bar, ...
                                              u, ...
                                              v_underline, ...
                                              v_bar, ... 
                                              x_spacing, ...
                                              y_spacing, ...
                                              bb_A_g, ...
                                              e_underline_g, ...
                                              e_bar_g);
    end

    %% Main Loop
    z_k_underline_trajectory(1, :) = [x_0_underline;
                                      d_0_underline]';
    z_k_bar_trajectory(1, :) = [x_0_bar;
                                d_0_bar]';
    if verbose
        log("[z_k_underline, z_k_bar]:");
        log([z_k_underline_trajectory(1, :)', z_k_bar_trajectory(1, :)']);
    end
                                
    x_k_minus_one_underline = x_0_underline;
    x_k_minus_one_bar = x_0_bar;
    d_k_minus_one_underline = d_0_underline;
    d_underline(1, :) = d_0_underline;
    d_k_minus_one_bar = d_0_bar;
    d_bar(1, :) = d_0_bar;
    u_k_minus_one = u(1, :)';
    z_k_minus_one_underline = [x_k_minus_one_underline;
                               d_k_minus_one_underline];
    z_k_minus_one_bar = [x_k_minus_one_bar;
                         d_k_minus_one_bar];
    xi_k_minus_one_underline = [z_k_minus_one_underline;
                                u_k_minus_one;
                                w_underline];
    xi_k_minus_one_bar = [z_k_minus_one_bar;
                          u_k_minus_one;
                          w_bar];
    h_k_minus_one_underline = h_0_underline;
    h_k_minus_one_bar = h_0_bar;
    z_k_minus_one_p_underline = [x_0_underline;
                                 d_0_underline];
    z_k_minus_one_p_bar = [x_0_bar;
                           d_0_bar];
    for k = 1:K
        if verbose
            log(sprintf("k: %u", k));
        end
        [bb_A_k_f, bb_W_k_f, bb_B_k_f, e_k_tilde_f] = ...
            get_observer_gains_f(f, xi_k_minus_one_underline, ...
                                 xi_k_minus_one_bar, ...
                                 num_grid_points_per_dim_f_domain, ...
                                 L_f, n, p, m, n_w, ...
                                 bb_A_f, e_bar_f, e_underline_f, k);
        [bb_A_k_h, bb_W_k_h, bb_B_k_h, e_k_tilde_h] = ...
            get_observer_gains_h(h_k_minus_one_underline, ...
                                 h_k_minus_one_bar, ...
                                 xi_k_minus_one_underline, ...
                                 xi_k_minus_one_bar, ...
                                 num_grid_points_per_dim_h_domain, ...
                                 L_h, n, p, m, n_w, ...
                                 bb_A_h, e_bar_h, e_underline_h, k);
        [x_k_a_p_underline, x_k_a_p_bar] = ...
            seven(bb_A_k_f, bb_B_k_f, bb_W_k_f, e_k_tilde_f, ...
                  z_k_minus_one_p_underline, z_k_minus_one_p_bar, ...
                  u_k_minus_one, w_underline, w_bar, n);
        [x_k_p_underline, x_k_p_bar] = ...
            four_a(f_d, z_k_minus_one_underline, z_k_minus_one_bar, ...
                   u_k_minus_one, w_underline, w_bar, ...
                   x_k_a_p_underline, x_k_a_p_bar);
        [d_k_p_underline, d_k_p_bar] = ...
            four_b(bb_A_k_h, bb_B_k_h, bb_W_k_h, e_k_tilde_h, ...
                   z_k_minus_one_p_underline, z_k_minus_one_p_bar, ...
                   u_k_minus_one, w_underline, w_bar, p);
        [z_k_p_underline, z_k_p_bar] = ...
            four_c(x_k_p_underline, x_k_p_bar, d_k_p_underline, d_k_p_bar);
        
        % Note: The following two lines are not part of the paper.
        z_k_p_underline = max([z_k_p_underline, z_0_underline], [], 2);
        z_k_p_bar = min([z_k_p_bar, z_0_bar], [], 2);
        
        [z_0_k_u_underline, z_0_k_u_bar] = eight(z_k_p_underline, ...
                                                  z_k_p_bar);
        z_i_minus_one_k_u_underline = z_0_k_u_underline;
        z_i_minus_one_k_u_bar = z_0_k_u_bar;
        y_k = y(k + 1, :)';
        u_k = u(k + 1, :)';
        for i = 1:num_interval_optimizations
            [A_i_k_g, B_i_k_g, W_i_k_g, e_i_k_underline_g, ...
             e_i_k_bar_g] = ...
                get_observer_gains_g(g, z_i_minus_one_k_u_underline, ...
                                     z_i_minus_one_k_u_bar, ...
                                     u_k_minus_one, v_underline, ...
                                     v_bar, ...
                                     num_grid_points_per_dim_g_domain, ...
                                     L_g, n, p, m, l, ...
                                     bb_A_g, e_bar_g, e_underline_g);
            [t_i_k_underline, t_i_k_bar] = ten(y_k, u_k, v_underline, ...
                                               v_bar, B_i_k_g, W_i_k_g, ...
                                               e_i_k_underline_g, ...
                                               e_i_k_bar_g, l);
            [alpha_i_k_underline, alpha_i_k_bar] = ...
                eleven(t_i_k_underline, t_i_k_bar, ...
                       z_i_minus_one_k_u_underline, ...
                       z_i_minus_one_k_u_bar, A_i_k_g);
            A_i_k_g_pinv = pinv(A_i_k_g);
            omega_i_k = calculate_omega(A_i_k_g, A_i_k_g_pinv, n, p);
            [z_i_k_u_underline, z_i_k_u_bar] = ...
                nine(A_i_k_g_pinv, alpha_i_k_underline, alpha_i_k_bar, ...
                omega_i_k, z_i_minus_one_k_u_underline, ...
                z_i_minus_one_k_u_bar);
            z_i_minus_one_k_u_underline = z_i_k_u_underline;
            z_i_minus_one_k_u_bar = z_i_k_u_bar;
        end
        [z_k_underline, z_k_bar] = five_a(z_i_k_u_underline, z_i_k_u_bar);
        [x_k_underline, x_k_bar, d_k_underline, d_k_bar] = ...
            five_b(z_k_underline, z_k_bar, n, p);
        d_underline(k + 1) = d_k_underline;
        d_bar(k + 1) = d_k_bar;
        xi_k_underline = [z_k_underline;
                          u_k;
                          w_underline];
        xi_k_bar = [z_k_bar;
                    u_k;
                    w_bar];
        xi_tilde(k + 1, :) = (1/2) * (xi_k_bar + xi_k_underline);
        epsilon(k + 1, :) = 2 * L_h * norm(xi_k_bar - xi_k_underline);
        h_k_bar = six_a(k, d_bar, L_h, xi_tilde, epsilon);
        h_k_underline = six_b(k, d_underline, L_h, xi_tilde, ...
                                      epsilon);
        z_k_underline_trajectory(k + 1, :) = [x_k_underline;
                                              d_k_underline]';
        z_k_bar_trajectory(k + 1, :) = [x_k_bar;
                                        d_k_bar]';
        if verbose
            log("[z_k_underline, z_k_bar]:");
            log([z_k_underline_trajectory(1, :)', ...
                 z_k_bar_trajectory(1, :)']);
        end
        
        % Prepare variables for next iteration.
        x_k_minus_one_underline = x_k_underline;
        x_k_minus_one_bar = x_k_bar;
        d_k_minus_one_underline = d_k_underline;
        d_k_minus_one_bar = d_k_bar;
        u_k_minus_one = u_k;
        z_k_minus_one_underline = [x_k_minus_one_underline;
                                   d_k_minus_one_underline];
        z_k_minus_one_bar = [x_k_minus_one_bar;
                             d_k_minus_one_bar];
        xi_k_minus_one_underline = xi_k_underline;
        xi_k_minus_one_bar = xi_k_bar;
        h_k_minus_one_underline = h_k_underline;
        h_k_minus_one_bar = h_k_bar;
        z_k_minus_one_p_underline = z_k_p_underline;
        z_k_minus_one_p_bar = z_k_p_bar;
    end
end


function [bb_A_k_f, bb_W_k_f, bb_B_k_f, e_k_tilde_f] = ...
    get_observer_gains_f(f, xi_k_minus_one_underline, ...
                         xi_k_minus_one_bar, ...
                         num_grid_points_per_dim_f_domain, ...
                         L_f, n, p, m, n_w, bb_A_f, ...
                         e_bar_f, e_underline_f, k)
    % Get the observer gains for f on the current interval.
    if k == 1
        % Case X = B in parallel affine abstraction.
        A_B_f = bb_A_f;
        e_k_bar_f = e_bar_f;
        e_k_underline_f = e_underline_f;
    else
        sigma_upper_bound_f = ...
            compute_sigma_upper_bound(L_f, xi_k_minus_one_underline, ...
                                      xi_k_minus_one_bar);
        [A_B_f, e_k_bar_f, e_k_underline_f] = ...
            local_parallel_affine_abstraction( ...
                f, f, xi_k_minus_one_underline, xi_k_minus_one_bar, ...
                num_grid_points_per_dim_f_domain, ...
                sigma_upper_bound_f, n, ...
                bb_A_f, e_bar_f, e_underline_f);
    end
    
    A_k_f = A_B_f(:, 1:n + p);
    B_k_f = A_B_f(:, n + p + 1:n + p + m);
    W_k_f = A_B_f(:, n + p + m + 1:n + p + m + n_w);
    
    A_k_f_plus = upper_plus(A_k_f);
    minus_A_k_f_plus_plus = - upper_plus_plus(A_k_f);
    bb_A_k_f = [A_k_f_plus, minus_A_k_f_plus_plus;
                minus_A_k_f_plus_plus, A_k_f_plus];
    W_k_f_plus = upper_plus(W_k_f);
    minus_W_k_f_plus_plus = - upper_plus_plus(W_k_f);
    bb_W_k_f = [W_k_f_plus, minus_W_k_f_plus_plus;
                minus_W_k_f_plus_plus, W_k_f_plus];
    bb_B_k_f = [B_k_f;
                B_k_f];
    e_k_tilde_f = [e_k_bar_f;
                   e_k_underline_f];
end


function [bb_A_k_h, bb_W_k_h, bb_B_k_h, e_k_tilde_h] = ...
    get_observer_gains_h(h_k_minus_one_underline, h_k_minus_one_bar, ...
                         xi_k_minus_one_underline, xi_k_minus_one_bar, ...
                         num_grid_points_per_dim_h_domain, ...
                         L_h, n, p, m, n_w, bb_A_h, ...
                         e_bar_h, e_underline_h, k)
    % Get the observer gains for h on the current interval.
    if k == 1
        % Case X = B in parallel affine abstraction.
        A_B_h = bb_A_h;
        e_k_bar_h = e_bar_h;
        e_k_underline_h = e_underline_h;
    else
        sigma_upper_bound_h = ...
            compute_sigma_upper_bound(L_h, xi_k_minus_one_underline, ...
                                      xi_k_minus_one_bar);
        [A_B_h, e_k_bar_h, e_k_underline_h] = ...
            local_parallel_affine_abstraction( ...
                h_k_minus_one_underline, h_k_minus_one_bar, ...
                xi_k_minus_one_underline, xi_k_minus_one_bar, ...
                num_grid_points_per_dim_h_domain, sigma_upper_bound_h, ...
                p, bb_A_h, e_bar_h, e_underline_h);
    end
    
    A_k_h = A_B_h(:, 1:n + p);
    B_k_h = A_B_h(:, n + p + 1:n + p + m);
    W_k_h = A_B_h(:, n + p + m + 1:n + p + m + n_w);
    
    A_k_h_plus = upper_plus(A_k_h);
    minus_A_k_h_plus_plus = - upper_plus_plus(A_k_h);
    bb_A_k_h = [A_k_h_plus, minus_A_k_h_plus_plus;
                minus_A_k_h_plus_plus, A_k_h_plus];
    W_k_h_plus = upper_plus(W_k_h);
    minus_W_k_h_plus_plus = - upper_plus_plus(W_k_h);
    bb_W_k_h = [W_k_h_plus, minus_W_k_h_plus_plus;
                minus_W_k_h_plus_plus, W_k_h_plus];
    bb_B_k_h = [B_k_h;
                B_k_h];
    e_k_tilde_h = [e_k_bar_h;
                   e_k_underline_h];
end


function [A_i_k_g, B_i_k_g, W_i_k_g, e_i_k_underline_g, e_i_k_bar_g] = ...
    get_observer_gains_g(g, z_i_minus_one_k_u_underline, ...
                         z_i_minus_one_k_u_bar, u_k_minus_one, ...
                         v_underline, v_bar, ...
                         num_grid_points_per_dim_g_domain, ...
                         L_g, n, p, m, l, bb_A_g, ...
                         e_bar_g, e_underline_g)
    % Get the observer gains for g on the current interval.
    xi_underline_g = [z_i_minus_one_k_u_underline;
                      u_k_minus_one;
                      v_underline];
    xi_bar_g = [z_i_minus_one_k_u_bar;
                u_k_minus_one;
                v_bar];
    sigma_upper_bound_g = compute_sigma_upper_bound(L_g, ...
                                                    xi_underline_g, ...
                                                    xi_bar_g);
    [A_B_g, e_i_k_bar_g, e_i_k_underline_g] = ...
        local_parallel_affine_abstraction( ...
            g, g, xi_underline_g, xi_bar_g, ...
            num_grid_points_per_dim_g_domain, sigma_upper_bound_g, l, ...
            bb_A_g, e_bar_g, e_underline_g);
    
    A_i_k_g = A_B_g(:, 1:n + p);
    B_i_k_g = A_B_g(:, n + p + 1:n + p + m);
    W_i_k_g = A_B_g(:, n + p + m + 1:n + p + m + l);
end


function [x_k_p_underline, x_k_p_bar] = ...
    four_a(f_d, z_k_minus_one_underline, z_k_minus_one_bar, ...
           u_k_minus_one, w_underline, w_bar, x_k_a_p_underline, ...
           x_k_a_p_bar)
    % (4a) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    xi_x = [z_k_minus_one_underline;
            u_k_minus_one;
            w_underline];
    xi_y = [z_k_minus_one_bar;
            u_k_minus_one;
            w_bar];
    x_k_p_underline = max(f_d(xi_x, xi_y), x_k_a_p_underline);
    
    xi_x = [z_k_minus_one_bar;
            u_k_minus_one;
            w_bar];
    xi_y = [z_k_minus_one_underline;
            u_k_minus_one;
            w_underline];
    x_k_p_bar = min(f_d(xi_x, xi_y), x_k_a_p_bar);
end


function [d_k_p_underline, d_k_p_bar] = ...
    four_b(bb_A_k_h, bb_B_k_h, bb_W_k_h, e_k_tilde_h, ...
           z_k_minus_one_p_underline, z_k_minus_one_p_bar, ...
           u_k_minus_one, w_underline, w_bar, p)
    % (4b) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    %
    % TODO:
    %   Why do we use z_k_minus_1_p instead of z_k_minus_1 here?
    solution = ...
        bb_A_k_h * [z_k_minus_one_p_bar; z_k_minus_one_p_underline] + ...
        bb_B_k_h * u_k_minus_one + ...
        bb_W_k_h * [w_bar; w_underline] + ...
        e_k_tilde_h;
    d_k_p_bar = solution(1:p);
    d_k_p_underline = solution(p + 1:2 * p);
end


function [z_k_p_underline, z_k_p_bar] = ...
    four_c(x_k_p_underline, x_k_p_bar, d_k_p_underline, d_k_p_bar)
    % (4c) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    z_k_p_underline = [x_k_p_underline; d_k_p_underline];
    z_k_p_bar = [x_k_p_bar; d_k_p_bar];
end


function [z_k_underline, z_k_bar] = five_a(z_i_k_u_underline, z_i_k_u_bar)
	% (5a) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    z_k_underline = z_i_k_u_underline;
    z_k_bar = z_i_k_u_bar;
end


function [x_k_underline, x_k_bar, d_k_underline, d_k_bar] = ...
    five_b(z_k_underline, z_k_bar, n, p)
    % (5b) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
	x_k_underline = z_k_underline(1:n);
    x_k_bar = z_k_bar(1:n);
    d_k_underline = z_k_underline(n + 1:n + p);
    d_k_bar = z_k_bar(n + 1:n + p);
end


function h_k_bar_handle = six_a(k, d_bar, L_h, xi_tilde, epsilon)
    % (6a) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    %
    % h_k_bar is generated for each dimension j in {1, ..., p}
    % simultaniously.
    %
    % Note:
    %   - In contrast to (6a) I include epsilon within the maximum. This
    %     is solves an issue with the indices there and is intuitively
    %     right.
    %   - In contrast to (6a) we align the data, s.t. 
    %     d_bar_k_minus_t is paired with xi_tilde_k_minus_t_minus_1
    %     and epsilon_k_minus_t_minus_1.
    %   - In contrast to (6a) we include the initial upper bound of the
    %     unknown dynamic d within the minimum.
    p = size(L_h, 1);
    d_0_bar = d_bar(1, :)';
    function d_k_bar = h_k_bar(xi_k)
        d_k_bar = zeros(p, 1);
        for j = 1:p
            d_k_bar(j, 1) = ...
                min(d_bar(2:k + 1, j) + L_h(j) * ...
                    vecnorm(xi_k' - xi_tilde(1:k, :), 2, 2) + ...
                    epsilon(1:k, j));
            % Note: The following line is not part of the paper.
            d_k_bar(j, 1) = min([d_k_bar(j, 1), d_0_bar(j)]);
        end
    end
    h_k_bar_handle = @h_k_bar;
end


function h_k_underline_handle = six_b(k, d_underline, L_h, xi_tilde, ...
                                      epsilon)
    % (6b) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    %
    % h_k_underline is generated for each dimension j in {1, ..., p}
    % simultaniously.
    %
    % Note:
    %   - In contrast to (6b) we include epsilon within the maximum. This
    %     is solves an issue with the indices there and is intuitively
    %     right.
    %   - In contrast to (6b) we align the data, s.t. 
    %     d_underline_k_minus_t is paired with xi_tilde_k_minus_t_minus_1
    %     and epsilon_k_minus_t_minus_1.
    %   - In contrast to (6a) we include the initial lower bound of the
    %     unknown dynamic d within the maximum.
    %   - In contrast to (6b) we substract epsilon. This corresponds to
    %     the original paper "Data-Driven Model Invalidation for Unknown 
    %     Lipschitz Continuous Systems via Abstraction", theorem 1.
    p = size(L_h, 1);
    d_0_underline = d_underline(1, :)';
    function d_k_underline = h_k_underline(xi_k)
        d_k_underline = zeros(p, 1);
        for j = 1:p
            d_k_underline(j, 1) = ...
                max(d_underline(2:k + 1, j) + L_h(j) * ...
                    vecnorm(xi_k' - xi_tilde(1:k, :), 2, 2) - ...
                    epsilon(1:k, j));
            % Note: The following line is not part of the paper.
            d_k_underline(j, 1) = max([d_k_underline(j, 1), ...
                                       d_0_underline(j)]);
        end
    end
    h_k_underline_handle = @h_k_underline;
end


function [x_k_a_p_underline, x_k_a_p_bar] = ...
    seven(bb_A_k_f, bb_B_k_f, bb_W_k_f, e_k_tilde_f, ...
          z_k_minus_one_p_underline, z_k_minus_one_p_bar, ...
          u_k_minus_one, w_underline, w_bar, n)
    % (7) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    %
    % TODO:
    %   Why do we use z_k_minus_1_p instead of z_k_minus_1 here?
    solution = ...
        bb_A_k_f * [z_k_minus_one_p_bar; z_k_minus_one_p_underline] + ...
        bb_B_k_f * u_k_minus_one + ...
        bb_W_k_f * [w_bar; w_underline] + ...
        e_k_tilde_f;
    x_k_a_p_bar = solution(1:n);
    x_k_a_p_underline = solution(n + 1:2 * n);
end


function [z_0_k_u_underline, z_0_k_u_bar] = eight(z_k_p_underline, ...
                                                  z_k_p_bar)
    % (8) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    z_0_k_u_underline = z_k_p_underline;
    z_0_k_u_bar = z_k_p_bar;
end


function [z_i_k_u_underline, z_i_k_u_bar] = ...
    nine(A_i_k_g_pinv, alpha_i_k_underline, alpha_i_k_bar, omega_i_k, ...
         z_i_minus_one_k_u_underline, z_i_minus_one_k_u_bar)
    % (9) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    A_i_k_g_pinv_plus = upper_plus(A_i_k_g_pinv);
    A_i_k_g_pinv_plus_plus = upper_plus_plus(A_i_k_g_pinv);
    z_i_k_u_underline = ...
        max(A_i_k_g_pinv_plus * alpha_i_k_underline - ...
                A_i_k_g_pinv_plus_plus * alpha_i_k_bar - ...
                omega_i_k, ...
            z_i_minus_one_k_u_underline);
    z_i_k_u_bar = ...
        min(A_i_k_g_pinv_plus * alpha_i_k_bar - ...
                A_i_k_g_pinv_plus_plus * alpha_i_k_underline + ...
                omega_i_k, ...
            z_i_minus_one_k_u_bar);
end


function omega_i_k = calculate_omega(A_i_k_g, A_i_k_g_pinv, n, p)
    % Defined after (11) in from "Interval Observers for Simultaneous 
    % State and Model Estimation of Partially Known Nonlinear Systems".
    kappa = Inf;
    omega_i_k = kappa * rowsupp(eye(n + p) - A_i_k_g_pinv * A_i_k_g);
end


function [t_i_k_underline, t_i_k_bar] = ten(y_k, u_k, v_underline, ...
                                            v_bar, B_i_k_g, W_i_k_g, ...
                                            e_i_k_underline_g, ...
                                            e_i_k_bar_g, l)
    % (10) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    W_i_k_g_plus = upper_plus(W_i_k_g);
    W_i_k_g_plus_plus = upper_plus_plus(W_i_k_g);
    solution = [y_k - B_i_k_g * u_k; y_k - B_i_k_g * u_k] + ...
               [W_i_k_g_plus_plus, - W_i_k_g_plus; - W_i_k_g_plus, ...
               W_i_k_g_plus_plus] * [v_bar; v_underline] - ...
               [e_i_k_underline_g; e_i_k_bar_g];
    t_i_k_bar = solution(1:l);
    t_i_k_underline = solution(l + 1:2 * l);
end


function [alpha_i_k_underline, alpha_i_k_bar] = ...
    eleven(t_i_k_underline, t_i_k_bar, z_i_minus_one_k_u_underline, ...
           z_i_minus_one_k_u_bar, A_i_k_g)
    % (11) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
    A_i_k_g_plus = upper_plus(A_i_k_g);
    A_i_k_g_plus_plus = upper_plus_plus(A_i_k_g);
    alpha_i_k_underline = ...
        max(t_i_k_underline, ...
            A_i_k_g_plus * z_i_minus_one_k_u_underline - ...
                A_i_k_g_plus_plus * z_i_minus_one_k_u_bar);
    alpha_i_k_bar = ...
        min(t_i_k_bar, ...
            A_i_k_g_plus * z_i_minus_one_k_u_bar - ...
                A_i_k_g_plus_plus * z_i_minus_one_k_u_underline);
end


function [r] = rowsupp(matrix)
    % rowsupp() from (II. Preliminaries) in "Interval Observers for 
    % Simultaneous State and Model  Estimation of Partially Known 
    % Nonlinear Systems".
    sz = size(matrix);
    p = sz(1);
    q = sz(2);
    r = ones(p, 1);
    zero = zeros(1, q);
    for row = 1:p
        if matrix(row, :) == zero
            r(row) = 0;
        end
    end
end


function [matrix_plus] = upper_plus(matrix)
    % Calculates max(a_i_j, 0) elementwise on matrix where a_i_j are the
    % elements of matrix.
    %
    % Naming with respect to "Interval Observers for Simultaneous State 
    % and Model Estimation of Partially Known Nonlinear Systems".
    matrix_plus = arrayfun(@max_zero, matrix);
end


function [maximum] = max_zero(value)
    % max(value, 0).
    maximum = max([value; 0]);
end


function [matrix_plus_plus] = upper_plus_plus(matrix)
    % Calculates - min(a_i_j, 0) elementwise on matrix where a_i_j are the
    % elements of matrix.
    %
    % Naming with respect to "Interval Observers for Simultaneous State 
    % and Model Estimation of Partially Known Nonlinear Systems".
    matrix_plus_plus = - arrayfun(@min_zero, matrix);
end


function [minimum] = min_zero(value)
    % min(value, 0).
    minimum = min([value; 0]);
end


function [] = func_plot_global_affine_abstraction_f(f, n, p, m, n_w, ...
                                                    x_0_underline, ...
                                                    x_0_bar, ...
                                                    d_0_underline, ...
                                                    d_0_bar, ...
                                                    u, ...
                                                    w_underline, ...
                                                    w_bar, ... 
                                                    x_spacing, ...
                                                    y_spacing, ...
                                                    bb_A_f, ...
                                                    e_underline_f, ...
                                                    e_bar_f)
    % Plot the global affine abstraction of f for dimension 1 and 2 of
    % the f input as x- and y-axis and output dimension one of f as z-axis.
    num_x_axis_coordinates = ...
        floor((x_0_bar(1) - x_0_underline(1)) / x_spacing);
    num_y_axis_coordinates = ...
        floor((x_0_bar(2) - x_0_underline(2)) / y_spacing);
    num_coords = num_x_axis_coordinates * num_y_axis_coordinates;
    x_1_axis_coords = linspace(x_0_underline(1), x_0_bar(1), ...
                        num_x_axis_coordinates);
    x_2_axis_coords = linspace(x_0_underline(2), x_0_bar(2), ...
                        num_y_axis_coordinates);
    [x_1_coords, x_2_coords] = meshgrid(x_1_axis_coords, ...
                                        x_2_axis_coords);
    x_1_coords = reshape(x_1_coords, num_coords, 1);
    x_2_coords = reshape(x_2_coords, num_coords, 1);
    remaining_x_coords = ones(num_coords, n - 2) .* ...
                         (x_0_underline(3:end) + ...
                          (x_0_bar(3:end) - x_0_underline(3:end))./2)';
    d_coords = ones(num_coords, p) .* ...
               (d_0_underline + (d_0_bar - d_0_underline)./2)';
    u_coords = ones(num_coords, m) .* u(1, :)';
    w_coords = ones(num_coords, n_w) .* ...
               (w_underline + (w_bar - w_underline)./2)';
    xi_coords = [x_1_coords, x_2_coords, remaining_x_coords, ...
                 d_coords, u_coords, w_coords];
    xi_coords_cell = num2cell(xi_coords', 1);
    
    f_of_xi_coords = cellfun(f, xi_coords_cell, 'UniformOutput', false);
    f_of_xi_coords = cell2mat(f_of_xi_coords);
    f_of_xi_coords = f_of_xi_coords';
    f_1_coords = f_of_xi_coords(:, 1);
    
    f_abstraction_underline_coords = (bb_A_f * xi_coords' + ...
                                      e_underline_f)';
    f_1_abstraction_underline_coords = ...
        f_abstraction_underline_coords(:, 1);
    
    f_abstraction_bar_coords = (bb_A_f * xi_coords' + e_bar_f)';
    f_1_abstraction_bar_coords = f_abstraction_bar_coords(:, 1);
    
	figure
    plot3(x_1_coords, x_2_coords, f_1_coords, ...
          'DisplayName', 'f_1');
	hold on
	plot3(x_1_coords, x_2_coords, f_1_abstraction_underline_coords, ...
          'DisplayName', 'A^{f} \cdot \xi + e^f_{underline}');
	plot3(x_1_coords, x_2_coords, f_1_abstraction_bar_coords, ...
          'DisplayName', 'A^{f} \cdot \xi_f + e^f_{bar}');
	hold off
    title('Global Parallel Affine Abstraction of f');
    xlabel('x_1');
    ylabel('x_2');
    zlabel('f_1');
    legend;
end


function [] = func_plot_global_affine_abstraction_g(g, n, p, m, l, ...
                                                    x_0_underline, ...
                                                    x_0_bar, ...
                                                    d_0_underline, ...
                                                    d_0_bar, ...
                                                    u, ...
                                                    v_underline, ...
                                                    v_bar, ... 
                                                    x_spacing, ...
                                                    y_spacing, ...
                                                    bb_A_g, ...
                                                    e_underline_g, ...
                                                    e_bar_g)
    % Plot the global affine abstraction of g for dimension 1 and 2 of
    % the g input as x- and y-axis and output dimension one of g as z-axis.
    num_x_axis_coordinates = ...
        floor((x_0_bar(1) - x_0_underline(1)) / x_spacing);
    num_y_axis_coordinates = ...
        floor((x_0_bar(2) - x_0_underline(2)) / y_spacing);
    num_coords = num_x_axis_coordinates * num_y_axis_coordinates;
    x_1_axis_coords = linspace(x_0_underline(1), x_0_bar(1), ...
                        num_x_axis_coordinates);
    x_2_axis_coords = linspace(x_0_underline(2), x_0_bar(2), ...
                        num_y_axis_coordinates);
    [x_1_coords, x_2_coords] = meshgrid(x_1_axis_coords, ...
                                        x_2_axis_coords);
    x_1_coords = reshape(x_1_coords, num_coords, 1);
    x_2_coords = reshape(x_2_coords, num_coords, 1);
    remaining_x_coords = ones(num_coords, n - 2) .* ...
                         (x_0_underline(3:end) + ...
                          (x_0_bar(3:end) - x_0_underline(3:end))./2)';
    d_coords = ones(num_coords, p) .* ...
               (d_0_underline + (d_0_bar - d_0_underline)./2)';
    u_coords = ones(num_coords, m) .* u(1, :)';
    v_coords = ones(num_coords, l) .* ...
               (v_underline + (v_bar - v_underline)./2)';
    xi_g_coords = [x_1_coords, x_2_coords, remaining_x_coords, ...
                 d_coords, u_coords, v_coords];
    xi_g_coords_cell = num2cell(xi_g_coords', 1);
    
    g_of_xi_g_coords = cellfun(g, xi_g_coords_cell, 'UniformOutput', ...
                               false);
    g_of_xi_g_coords = cell2mat(g_of_xi_g_coords);
    g_of_xi_g_coords = g_of_xi_g_coords';
    g_1_coords = g_of_xi_g_coords(:, 1);
    
    g_abstraction_underline_coords = (bb_A_g * xi_g_coords' + ...
                                      e_underline_g)';
    g_1_abstraction_underline_coords = ...
        g_abstraction_underline_coords(:, 1);
    
    g_abstraction_bar_coords = (bb_A_g * xi_g_coords' + e_bar_g)';
    g_1_abstraction_bar_coords = g_abstraction_bar_coords(:, 1);
    
	figure
    plot3(x_1_coords, x_2_coords, g_1_coords, ...
          'DisplayName', 'g_1');
	hold on
	plot3(x_1_coords, x_2_coords, g_1_abstraction_underline_coords, ...
          'DisplayName', 'A^{g} \cdot \xi_g + e^g_{underline}');
	plot3(x_1_coords, x_2_coords, g_1_abstraction_bar_coords, ...
          'DisplayName', 'A^{g} \cdot \xi_g + e^g_{bar}');
	hold off
    title('Global Parallel Affine Abstraction of g');
    xlabel('x_1');
    ylabel('x_2');
    zlabel('g_1');
    legend;
end


function [] = func_plot_global_affine_abstraction_h(h, h_underline, ...
                                                    h_bar, ...
                                                    n, p, m, n_w, ...
                                                    x_0_underline, ...
                                                    x_0_bar, ...
                                                    d_0_underline, ...
                                                    d_0_bar, ...
                                                    u, ...
                                                    w_underline, ...
                                                    w_bar, ... 
                                                    x_spacing, ...
                                                    y_spacing, ...
                                                    bb_A_h, ...
                                                    e_underline_h, ...
                                                    e_bar_h)
    % Plot the global affine abstraction of h for dimension 1 and 2 of
    % the h input as x- and y-axis and output dimension one of h as z-axis.
    num_x_axis_coordinates = ...
        floor((x_0_bar(1) - x_0_underline(1)) / x_spacing);
    num_y_axis_coordinates = ...
        floor((x_0_bar(2) - x_0_underline(2)) / y_spacing);
    num_coords = num_x_axis_coordinates * num_y_axis_coordinates;
    x_1_axis_coords = linspace(x_0_underline(1), x_0_bar(1), ...
                        num_x_axis_coordinates);
    x_2_axis_coords = linspace(x_0_underline(2), x_0_bar(2), ...
                        num_y_axis_coordinates);
    [x_1_coords, x_2_coords] = meshgrid(x_1_axis_coords, ...
                                        x_2_axis_coords);
    x_1_coords = reshape(x_1_coords, num_coords, 1);
    x_2_coords = reshape(x_2_coords, num_coords, 1);
    remaining_x_coords = ones(num_coords, n - 2) .* ...
                         (x_0_underline(3:end) + ...
                          (x_0_bar(3:end) - x_0_underline(3:end))./2)';
    d_coords = ones(num_coords, p) .* ...
               (d_0_underline + (d_0_bar - d_0_underline)./2)';
    u_coords = ones(num_coords, m) .* u(1, :)';
    w_coords = ones(num_coords, n_w) .* ...
               (w_underline + (w_bar - w_underline)./2)';
    xi_coords = [x_1_coords, x_2_coords, remaining_x_coords, ...
                 d_coords, u_coords, w_coords];
    xi_coords_cell = num2cell(xi_coords', 1);
    
    if(isa(h, 'function_handle'))
        h_of_xi_coords = cellfun(h, xi_coords_cell, ...
                                 'UniformOutput', false);
        h_of_xi_coords = cell2mat(h_of_xi_coords);
        h_of_xi_coords = h_of_xi_coords';
        h_1_coords = h_of_xi_coords(:, 1);
    end
    
    h_underline_of_xi_coords = cellfun(h_underline, xi_coords_cell, ...
                             'UniformOutput', false);
    h_underline_of_xi_coords = cell2mat(h_underline_of_xi_coords);
    h_underline_of_xi_coords = h_underline_of_xi_coords';
    h_underline_1_coords = h_underline_of_xi_coords(:, 1);
    
    h_bar_of_xi_coords = cellfun(h_bar, xi_coords_cell, ...
                             'UniformOutput', false);
    h_bar_of_xi_coords = cell2mat(h_bar_of_xi_coords);
    h_bar_of_xi_coords = h_bar_of_xi_coords';
    h_bar_1_coords = h_bar_of_xi_coords(:, 1);
    
    h_abstraction_underline_coords = (bb_A_h * xi_coords' + ...
                                      e_underline_h)';
    h_1_abstraction_underline_coords = ...
        h_abstraction_underline_coords(:, 1);
    
    h_abstraction_bar_coords = (bb_A_h * xi_coords' + e_bar_h)';
    h_1_abstraction_bar_coords = h_abstraction_bar_coords(:, 1);
    
    figure
    plot3(x_1_coords, x_2_coords, h_underline_1_coords, ...
          'DisplayName', 'h_{underline, 1}');
	hold on
    if(isa(h, 'function_handle'))
        plot3(x_1_coords, x_2_coords, h_1_coords, ...
              'DisplayName', 'h_{1}');
    end
    plot3(x_1_coords, x_2_coords, h_bar_1_coords, ...
          'DisplayName', 'h_{bar, 1}');
	plot3(x_1_coords, x_2_coords, h_1_abstraction_underline_coords, ...
          'DisplayName', 'A^{h} \cdot \xi + e^h_{underline}');
	plot3(x_1_coords, x_2_coords, h_1_abstraction_bar_coords, ...
          'DisplayName', 'A^{h} \cdot \xi + e^h_{bar}');
	hold off
    title('Global Parallel Affine Abstraction of h');
    xlabel('x_1');
    ylabel('x_2');
    zlabel('d_1');
    legend;
end
    
