function [A, b, last_defined_A_row, last_defined_b_row] = ...
            add_twelve_b(xi, A, b, last_defined_A_row, ...
                         last_defined_b_row, blackboard_bold_A, e_bar, ...
                         e_underline, m_prime, n_prime)
    % (12b) from "Interval Observers for Simultaneous State and Model 
    % Estimation of Partially Known Nonlinear Systems".
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
