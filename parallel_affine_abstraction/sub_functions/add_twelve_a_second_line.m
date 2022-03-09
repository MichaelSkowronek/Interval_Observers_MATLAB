function [A, b, last_defined_A_row, last_defined_b_row] ...
    = add_twelve_a_second_line(A, b, last_defined_A_row, ...
                               last_defined_b_row, m_prime, sigma)
    % (12a) second line from "Interval Observers for Simultaneous State 
    % and Model Estimation of Partially Known Nonlinear Systems".
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
