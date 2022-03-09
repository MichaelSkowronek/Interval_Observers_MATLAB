function [A, b, last_defined_A_row, last_defined_b_row] = ...
    add_twelve_a_first_line(xi, A, b, last_defined_A_row, ...
                            last_defined_b_row, q_underline, q_bar, ...
                            m_prime, n_prime, sigma)
    % (12a) first line from "Interval Observers for Simultaneous State 
    % and Model Estimation of Partially Known Nonlinear Systems".
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
