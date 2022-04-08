function [A_B, e_B_bar, e_B_underline] = ...
    global_parallel_affine_abstraction( ...
    q_underline, q_bar, xi_underline, xi_bar, ...
    num_grid_points_per_dim, sigma_upper_bound, m_prime)
    % Proposition 2 (Parallel Affine Abstraction) from the paper "Interval 
    % Observers for Simultaneous State and Model Estimation of Partially 
    % Known Nonlinear Systems" to calculate blackboard
    % bold A, e_bar, e_underline for global affinity, i.e. B = X.
    %
    % Note: The power index 'q' is left out from variable naming.
    % Note: The 's' index in xi_s is left out.
    %
    % Args:
    %   sigma_upper_bound: Instead of sigma, we also can provide an 
    %                      upper bound. Size (m_prime, 1).
    %   m_prime: The output dimension of q_underline and q_bar.
    %
    % Returns:
    %   A_B: A_B = blackboard_bold_A in this case. Size (m_prime, n_prime).
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
