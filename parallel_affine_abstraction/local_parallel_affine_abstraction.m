function [A_B, e_B_bar, e_B_underline] = ...
    local_parallel_affine_abstraction( ...
    q_underline, q_bar, xi_underline, xi_bar, ...
    num_grid_points_per_dim, sigma_upper_bound, m_prime, ...
    blackboard_bold_A, e_bar, e_underline)
    % Proposition 2 (Parallel Affine Abstraction) from the paper "Interval 
    % Observers for Simultaneous State and Model Estimation of Partially 
    % Known Nonlinear Systems" to calculate A_B, e_B_bar, e_B_underline 
    % for local affinity.
    %
    % Note: The power index 'q' is left out from variable naming.
    % Note: The 's' index in xi_s is left out.
    %
    % Args:
    %   sigma_upper_bound: Instead of sigma we also can provide an upper 
    %                      bound. Size (m_prime, 1).
    %   m_prime: The output dimension of q_underline and q_bar.
    %   blackboard_bold_A: The global parallel affine abstraction matrix
    %                      A_B with B = X.
    %   e_bar: The global parallel affine abstraction vector e_B_bar with
    %          B = X.
    %   e_underline: The global parallel affine abstraction vector 
    %                e_B_underline with B = X.
    %
    % Returns:
    %   A_B: Size (m_prime, n_prime).
    %   e_B_bar: Size (m_prime, 1).
    %   e_B_underline: Size (m_prime, 1).
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
