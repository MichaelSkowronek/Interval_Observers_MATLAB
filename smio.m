function [] = ...
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
    %   f_d: Mixed-monotone mapping of f. Assume f: R^n -> R^m, then
    %        f_d: R^n x R^n -> R^m.
    %   y: A matrix of size (K + 1 x l) of obervations.
    %      The first entry is the observation of step 0.
    %   u: A matrix of size (K + 1 x m) of input vectors.
    %      The first entry is the input of step 0.
    %   L_f: The Lipschitz constant of f.
    %   L_g: The Lipschitz constant of g.
    %   L_h: A vector of size p of Lipschitz constants for each output
    %        dimension of h.
    %   K: The number of observations and corresponding iterations k.
    %   num_interval_optimizations: The number of optimizations i of the 
    %                               interval within each observation 
    %                               step k.
    % varargin:
    %   Key-Word Args:
    %       plot_global_affine_abstraction_f: Default false, type boolean.
    %       plot_global_affine_abstraction_g: Default false, type boolean.
    %       plot_global_affine_abstraction_h: Default false, type boolean.
    %       x_spacing: Default 0.01. The spacing between plotting grid 
    %                  points on the x-axis.
    %       y_spacing: Default 0.01. The spacing between plotting grid 
    %                  points on the y-axis.
    %
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
    addParameter(input_parser, 'x_spacing', 0.01);
    addParameter(input_parser, 'y_spacing', 0.01);
    parse(input_parser, varargin{:});
    results = input_parser.Results;
    plot_global_affine_abstraction_f = ...
        results.plot_global_affine_abstraction_f;
    plot_global_affine_abstraction_h = ...
        results.plot_global_affine_abstraction_h;
    plot_global_affine_abstraction_g = ...
        results.plot_global_affine_abstraction_g;
    x_spacing = results.x_spacing;
    y_spacing = results.y_spacing;
    
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
    lambda_h = max(L_h);  % TODO: Does this make sense?
    k = 0;
    h_0_bar = six_a(k, d_bar, L_h, xi_tilde, epsilon);
    h_0_underline = six_b(k, d_underline, L_h, xi_tilde, epsilon);
    sigma_upper_bound_h = compute_sigma_upper_bound(lambda_h, ...
                                                  xi_0_underline, ...
                                                  xi_0_bar);
    [bb_A_h, e_bar_h, e_underline_h] = ...
        global_parallel_affine_abstraction( ...
            h_0_underline, h_0_bar, xi_0_underline, xi_0_bar, ...
            num_grid_points_per_dim_h_domain, sigma_upper_bound_h, p);
        
    if plot_global_affine_abstraction_h
        func_plot_global_affine_abstraction_h(h_0_underline, h_0_bar, ...
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
    % TODO: Doing this here is not wrong, but is there a better way?
    %       We maybe could do it within each iteration.
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
    % TODO: Is this the right way to calculate z_0_p_underline and
    %       z_0_p_bar.
    z_k_minus_one_p_underline = [x_0_underline;
                                 d_0_underline];
    z_k_minus_one_p_bar = [x_0_bar;
                                 d_0_bar];
    for k = 1:K
        [bb_A_k_f, bb_W_k_f, bb_B_k_f, e_k_tilde_f] = ...
            get_observer_gains_f(f, xi_k_minus_one_underline, ...
                                 xi_k_minus_one_bar, ...
                                 num_grid_points_per_dim_f_domain, ...
                                 sigma_upper_bound_f, n, p, m, n_w, ...
                                 bb_A_f, e_bar_f, e_underline_f, k);
        [bb_A_k_h, bb_W_k_h, bb_B_k_h, e_k_tilde_h] = ...
            get_observer_gains_h(h_k_minus_one_underline, ...
                                 h_k_minus_one_bar, ...
                                 xi_k_minus_one_underline, ...
                                 xi_k_minus_one_bar, ...
                                 num_grid_points_per_dim_h_domain, ...
                                 sigma_upper_bound_h, n, p, m, n_w, ...
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
                                     sigma_upper_bound_g, n, p, m, l, ...
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
            omega_i_k = calculate_omega(kappa, A_i_k_g, A_i_k_g_pinv, ...
                                        n, p);
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
                         sigma_upper_bound_f, n, p, m, n_w, bb_A_f, ...
                         e_bar_f, e_underline_f, k)
    % Get the observer gains for f on the current interval.
    if k == 1
        % Case X = B in parallel affine abstraction.
        A_B_f = bb_A_f;
        e_k_bar_f = e_bar_f;
        e_k_underline_f = e_underline_f;
    else
        [A_B_f, e_k_bar_f, e_k_underline_f] = ...
            local_parallel_affine_abstraction( ...
                f, f, xi_k_minus_one_underline, xi_k_minus_one_bar, ...
                num_grid_points_per_dim_f_domain, sigma_upper_bound_f, n, ...
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
                         sigma_upper_bound_h, n, p, m, n_w, bb_A_h, ...
                         e_bar_h, e_underline_h, k)
    % Get the observer gains for h on the current interval.
    if k == 1
        % Case X = B in parallel affine abstraction.
        A_B_h = bb_A_h;
        e_k_bar_h = e_bar_h;
        e_k_underline_h = e_underline_h;
    else
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
                         sigma_upper_bound_g, n, p, m, l, bb_A_g, ...
                         e_bar_g, e_underline_g)
    % Get the observer gains for g on the current interval.
    xi_underline_g = [z_i_minus_one_k_u_underline;
                      u_k_minus_one;
                      v_underline];
    xi_bar_g = [z_i_minus_one_k_u_bar;
                u_k_minus_one;
                v_bar];
    [A_B_g, e_i_k_bar_g, e_i_k_underline_g] = ...
        local_parallel_affine_abstraction( ...
            g, g, xi_underline_g, xi_bar_g, ...
            num_grid_points_per_dim_g_domain, sigma_upper_bound_g, l, ...
            bb_A_g, e_bar_g, e_underline_g);
    
    A_i_k_g = A_B_g(:, 1:n + p);
    B_i_k_g = A_B_g(:, n + p + 1:n + p + m);
    W_i_k_g = A_B_g(:, n + p + m + 1:n + p + m + l);
end

    
function [A_B, e_B_bar, e_B_underline] = ...
    global_parallel_affine_abstraction( ...
    q_underline, q_bar, xi_underline, xi_bar, ...
    num_grid_points_per_dim, sigma_upper_bound, m_prime)
    % Proposition 2 (Parallel Affine Abstraction) to calculate blackboard
    % bold A, e_bar, e_underline for global affinity, i.e. B = X.
    %
    % Note: The power index 'q' is left out from variable naming.
    % Note: The 's' index in xi_s is left out.
    %
    % Args:
    %   sigma_upper_bound: Instead of sigma we also can give an upper 
    %                      bound.
    %   m_prime: The output dimension of q_underline and q_bar.
    %
    % Returns:
    %   A_B: A_B = bb(A) in this case. Size (m_prime, n_prime).
    %   e_B_bar: e_B_bar = e_bar in this case. Size (m_prime, 1).
    %   e_B_underline: e_B_underline = e_underline in this case. Size 
    %                  (m_prime, 1).
    sigma = sigma_upper_bound;
    n_prime = size(xi_underline, 1);
    V_B = generate_V_B(xi_underline, xi_bar, num_grid_points_per_dim);
    num_grid_points = size(V_B, 1);
    % Note: Notation follows the definition of linprog() from now on.
    b = zeros(m_prime + num_grid_points * m_prime, 1);
    last_defined_b_row = 0;
    last_defined_A_row = 0;
    
    % Define x.
    % x = concat(theta, flatten(A_B), e_B_bar, e_B_underline).
    A_B_num_params = m_prime * n_prime;
    x_length = 1 + A_B_num_params + 2 * m_prime;
    A = zeros(m_prime + num_grid_points * m_prime, x_length);
    
    % Define f.
    f = zeros(x_length, 1);
    f(1) = 1;  % We only want theta to be minimized.
    
    % (12a) second line.
    [A, b, last_defined_A_row, last_defined_b_row] = ...
        add_twelve_a_second_line(A, b, last_defined_A_row, ...
                                 last_defined_b_row, m_prime, sigma);
    
    % (12a) first line.
    for i = 1:num_grid_points
        xi = V_B(i, :)';
        [A, b, last_defined_A_row, last_defined_b_row] = ...
            add_twelve_a_first_line(xi, A, b, last_defined_A_row, ...
                                    last_defined_b_row, q_underline, ...
                                    q_bar, m_prime, n_prime, sigma);
    end
    
    x = linprog(f, A, b);
    A_B = reshape(x(1 + 1:1 + m_prime * n_prime), n_prime, m_prime)';
    e_B_bar = x(end - 2 * m_prime + 1:end - m_prime);
    e_B_underline = x(end - m_prime + 1:end);
end


function [A_B, e_B_bar, e_B_underline] = ...
    local_parallel_affine_abstraction( ...
    q_underline, q_bar, xi_underline, xi_bar, ...
    num_grid_points_per_dim, sigma_upper_bound, m_prime, ...
    blackboard_bold_A, e_bar, e_underline)
    % Proposition 2 (Parallel Affine Abstraction) to calculate blackboard
    % bold A, e_bar, e_underline for global affinity, i.e. B = X.
    %
    % Note: The power index 'q' is left out from variable naming.
    % Note: The 's' index in xi_s is left out.
    %
    % Args:
    %   sigma_upper_bound: Instead of sigma we also can give an upper 
    %                      bound.
    %   m_prime: The output dimension of q_underline and q_bar.
    %   blackboard_bold_A: The global parallel affine abstraction matrix
    %                      A_B with B = X.
    %   e_bar: The global parallel affine abstraction vector e_B_bar with
    %          B = X.
    %   e_underline: The global parallel affine abstraction vector 
    %                e_B_underline with B = X.
    %
    % Returns:
    %   A_B: A_B = bb(A) in this case. Size (m_prime, n_prime).
    %   e_B_bar: e_B_bar = e_bar in this case. Size (m_prime, 1).
    %   e_B_underline: e_B_underline = e_underline in this case. Size 
    %                  (m_prime, 1).
    sigma = sigma_upper_bound;
    n_prime = size(xi_underline, 1);
    V_B = generate_V_B(xi_underline, xi_bar, num_grid_points_per_dim);
    num_grid_points = size(V_B, 1);
    % Note: Notation follows the definition of linprog() from now on.
    b = zeros(m_prime + num_grid_points * 2 * m_prime, 1);
    last_defined_b_row = 0;
    last_defined_A_row = 0;
    
    % Define x.
    % x = concat(theta, flatten(A_B), e_B_bar, e_B_underline).
    A_B_num_params = m_prime * n_prime;
    x_length = 1 + A_B_num_params + 2 * m_prime;
    A = zeros(m_prime + num_grid_points * 2 * m_prime, x_length);
    
    % Define f.
    f = zeros(x_length, 1);
    f(1) = 1;  % We only want theta to be minimized.
    
    % (12a) second line.
    [A, b, last_defined_A_row, last_defined_b_row] = ...
        add_twelve_a_second_line(A, b, last_defined_A_row, ...
                                 last_defined_b_row, m_prime, sigma);
    
    % (12a) first line and (12b).
    for i = 1:num_grid_points
        xi = V_B(i, :)';
        [A, b, last_defined_A_row, last_defined_b_row] = ...
            add_twelve_a_first_line(xi, A, b, last_defined_A_row, ...
                                    last_defined_b_row, q_underline, ...
                                    q_bar, m_prime, n_prime, sigma);
        [A, b, last_defined_A_row, last_defined_b_row] = ...
            add_twelve_b(xi, A, b, last_defined_A_row, ...
                         last_defined_b_row, blackboard_bold_A, e_bar, ...
                         e_underline, m_prime, n_prime);
    end
    
    x = linprog(f, A, b);
    A_B = reshape(x(1 + 1:1 + m_prime * n_prime), n_prime, m_prime)';
    e_B_bar = x(end - 2 * m_prime + 1:end - m_prime);
    e_B_underline = x(end - m_prime + 1:end);
end


function [A, b, last_defined_A_row, last_defined_b_row] ...
    = add_twelve_a_second_line(A, b, last_defined_A_row, ...
                               last_defined_b_row, m_prime, sigma)
    % (12a) second line.
    % Add e_B_bar - e_B_underline - 2*sigma <= theta * 1_{m_prime}
    % <=> 1 * e_B_bar + (-1) * e_B_underline + (-1) * theta * 1_{m_prime} 
    %     <= 2*sigma.
    % Add (2*sigma).
    b(last_defined_b_row + 1:last_defined_b_row + m_prime) = 2 * sigma;
    last_defined_b_row = last_defined_b_row + m_prime;
    for i = 1:m_prime
        % Add (+ e_B_bar).
        A(last_defined_A_row + 1, end - 2 * m_prime + i) = 1;
        % Add (+ (-1) * e_B_underline).
        A(last_defined_A_row + 1, end - m_prime + i) = -1;
        % Add (+ (-1) * theta * 1_{m_prime}).
        A(last_defined_A_row + 1, 1) = -1;
        last_defined_A_row = last_defined_A_row + 1;
    end
end


function [A, b, last_defined_A_row, last_defined_b_row] = ...
    add_twelve_a_first_line(xi, A, b, last_defined_A_row, ...
                            last_defined_b_row, q_underline, q_bar, ...
                            m_prime, n_prime, sigma)
    % (12a) first line.
    q_underline_xi = q_underline(xi);
    q_bar_xi = q_bar(xi);

    % Add A_B * xi + e_B_underline + sigma <= q_underline(xi)
    % <=> A_B * xi + e_B_underline <= 
    %     q_underline(xi) - sigma for each xi in V_B.
    % Add (q_underline(xi) - sigma).
    b(last_defined_b_row + 1: ...
      last_defined_b_row + m_prime) = q_underline_xi - sigma;
    last_defined_b_row = last_defined_b_row + m_prime;

    for j = 1:m_prime
        % Add (+ e_B_underline).
        A(last_defined_A_row + 1, end - m_prime + j) = 1;
        % Add (+ A_B * xi).
        A(last_defined_A_row + 1, ...
          1 + (j - 1) * n_prime + 1: ...
          1 + (j - 1) * n_prime + n_prime) = xi';
        last_defined_A_row = last_defined_A_row + 1;
    end

    % Add q_bar(xi) <= A_B * xi + e_B_bar - sigma
    % <=> - A_B * xi - e_B_bar <= - sigma - q_bar(xi).
    % Add (- sigma - q_bar(xi)).
    b(last_defined_b_row + 1: ...
      last_defined_b_row + m_prime) = - q_bar_xi - sigma;
    last_defined_b_row = last_defined_b_row + m_prime;

    for j = 1:m_prime
        % Add (- e_B_bar).
        A(last_defined_A_row + 1, end - 2 * m_prime + j) = - 1;
        % Add (- A_B * xi).
        A(last_defined_A_row + 1, ...
          1 + (j - 1) * n_prime + 1: ...
          1 + (j - 1) * n_prime + n_prime) = - xi';
        last_defined_A_row = last_defined_A_row + 1;
    end
end


function [A, b, last_defined_A_row, last_defined_b_row] = ...
            add_twelve_b(xi, A, b, last_defined_A_row, ...
                         last_defined_b_row, blackboard_bold_A, e_bar, ...
                         e_underline, m_prime, n_prime)
    % (12b).
    bb_A_times_xi = blackboard_bold_A * xi;
    % Add e_underline - e_B_underline <= (A_B - blackboard_bold_A) * xi
    % <=> - A_B * xi - e_B_underline <= - e_underline - blackboard_bold_A
    %                                                 * xi.
    % Add (- e_underline - blackboard_bold_A * xi).
    b(last_defined_b_row + 1: ...
      last_defined_b_row + m_prime) = - e_underline - bb_A_times_xi;
    last_defined_b_row = last_defined_b_row + m_prime;
    
    for j = 1:m_prime
        % Add (- e_B_underline).
        A(last_defined_A_row + 1, end - m_prime + j) = - 1;
        % Add (- A_B * xi).
        A(last_defined_A_row + 1, ...
          1 + (j - 1) * n_prime + 1: ...
          1 + (j - 1) * n_prime + n_prime) = - xi';
        last_defined_A_row = last_defined_A_row + 1;
    end
    
    % Add (A_B - blackboard_bold_A) * xi <= e_bar - e_B_bar
    % <=> A_B * xi + e_B_bar <= e_bar + blackboard_bold_A * xi.
    % Add (e_bar + blackboard_bold_A * xi).
    b(last_defined_b_row + 1: ...
      last_defined_b_row + m_prime) = e_bar + bb_A_times_xi;
    last_defined_b_row = last_defined_b_row + m_prime;
    
    for j = 1:m_prime
        % Add (+ e_B_bar).
        A(last_defined_A_row + 1, end - 2 * m_prime + j) = 1;
        % Add (+ A_B * xi).
        A(last_defined_A_row + 1, ...
          1 + (j - 1) * n_prime + 1: ...
          1 + (j - 1) * n_prime + n_prime) = xi';
        last_defined_A_row = last_defined_A_row + 1;
    end
end


function [sigma_upper_bound] = compute_sigma_upper_bound(lambda, ...
                                                         xi_underline, ...
                                                         xi_bar)
    % Paper "Mesh-Based Affine Abstraction of Nonlinear Systems with 
    % Tighter Bounds", Proposition 1, case (ii).
    % 
    % Args:
    %    lambda: Lipschitz constant of the function f, see paper.
    %    xi_underline: The minimal value a verticy of S should have.
    %                  Shape (dimension, 1).
    %    xi_bar: The maximal value a verticy of S should have.
    %            Shape (dimension, 1). 
    % 
    % Returns:
    %    sigma_upper_bound: The case (ii) upper bound of sigma.
    % TODO:
    %   - Check if this function holds for zero width intervals still.
    n_plus_m = size(xi_underline, 1);
    
    % Count number of intervals with zero width
    num_zero_intervals = 0;
    width = xi_bar - xi_underline;
    for i = 1:n_plus_m
        if width(i) == 0
            num_zero_intervals = num_zero_intervals + 1;
        end
    end
    
    % Choose affine independent verticies of S.
    num_verticies = n_plus_m + 1 - num_zero_intervals;
    S_verticies = zeros(num_verticies, n_plus_m);
    S_verticies_offset_unit_vectors = zeros(num_verticies, n_plus_m);
    for i = 1:num_verticies
        is_affine_independent = false;
        while is_affine_independent == false
            vertex_i = rand(n_plus_m, 1);
            vertex_i = xi_underline + (xi_bar - xi_underline) .* vertex_i;
            
            % Check affine independence.
            if i == 1
                is_affine_independent = true;
            else
                vertex_i_offset = vertex_i - S_verticies(1, :)';
                vertex_i_offset_unit_vector = ...
                    vertex_i_offset ./ norm(vertex_i_offset);
                S_verticies_offset_unit_vectors(i, :) = ...
                    vertex_i_offset_unit_vector';
                matrix_rank = ...
                    rank(S_verticies_offset_unit_vectors(2:i, :));
                if matrix_rank == i - 1
                   is_affine_independent = true;
                end
            end
        end
        S_verticies(i, :) = vertex_i';
    end
    
    % Calculate diameter delta.
    delta = 0;
    for i = 1:num_verticies
       for j =  i+1:num_verticies
           distance = norm(S_verticies(i, :) - S_verticies(j, :));
           if distance > delta
              delta = distance;
           end
       end
    end
    
    % Calculate delta_s.
    delta_s = sqrt(n_plus_m/(2*(n_plus_m + 1)))*delta;
    
    % Calculate upper bound.
    sigma_upper_bound = lambda * delta_s;
end


function [V_B] = generate_V_B(xi_underline, xi_bar, ...
                              num_grid_points_per_dim)
    % Generates a regular grid of the hyperrectangle defined by 
    % xi_underline and xi_bar.
    %
    % Returns:
    %   V_B: Shape (num_grid_points, xi_dimensionality)
    n_prime = size(xi_underline, 1);
    
    % Calculate V_B vectors.
    num_grid_vectors = num_grid_points_per_dim^n_prime;
    V_B = zeros(num_grid_vectors, n_prime);
    for dimension = 1:n_prime
        smalles_value = xi_underline(dimension, 1);
        largest_value = xi_bar(dimension, 1);
        values_of_this_dim = linspace(smalles_value, largest_value, ...
                                      num_grid_points_per_dim)';
        num_equal_values_in_a_row = num_grid_points_per_dim^(n_prime - ...
                                                             dimension);
        intermediate_matrix = repmat(values_of_this_dim, 1, ...
                                     num_equal_values_in_a_row);
        equal_values_in_a_row = reshape(intermediate_matrix', [], 1);
        num_repeats_of_vector = num_grid_points_per_dim^(dimension - 1);
        V_B(:, dimension) = repmat(equal_values_in_a_row, ...
                                   num_repeats_of_vector, 1);
    end
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
    %   The paper does not include epsilon within the parenthesis as I do
    %   here. But without it makes no sense.
    p = size(L_h, 1);
    function d_k_bar = h_k_bar(xi_k)
        d_k_bar = zeros(p, 1);
        for j = 1:p
            d_k_bar(j, 1) = min(d_bar(1:k + 1, j) + ...
                            L_h(j) * ...
                            norm(xi_k' - xi_tilde(1:k + 1, :)) + ...
                            epsilon(1:k + 1, j));
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
    %   The paper does not include epsilon within the parenthesis as I do
    %   here. But without it makes no sense.
    p = size(L_h, 1);
    function d_k_underline = h_k_underline(xi_k)
        d_k_underline = zeros(p, 1);
        for j = 1:p
            d_k_underline(j, 1) = max(d_underline(1:k + 1, j) + ...
                                  L_h(j) * ...
                                  norm(xi_k' - xi_tilde(1:k + 1, :)) + ...
                                  epsilon(1:k + 1, j));
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
    %
    % TODO:
    %   - Test this function.
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
    
    plot3(x_1_coords, x_2_coords, f_1_coords, ...
          x_1_coords, x_2_coords, f_1_abstraction_underline_coords, ...
          x_1_coords, x_2_coords, f_1_abstraction_bar_coords);
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
    zlabel('d_1');
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
    zlabel('d_1');
    legend;
end


function [] = func_plot_global_affine_abstraction_h(h_underline, ...
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
    
