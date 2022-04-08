function [sigma_upper_bound] = compute_sigma_upper_bound(lambda, ...
                                                         xi_underline, ...
                                                         xi_bar)
    % Paper "Mesh-Based Affine Abstraction of Nonlinear Systems with 
    % Tighter Bounds", Proposition 1, case (ii).
    %
    % If the corresponding function f is a vector field, a
    % sigma_upper_bound for each codomain dimension can be computet by
    % providing a lambda vector with the codomain dimension. These are
    % independent of each other.
    % 
    % Args:
    %    lambda: Lipschitz constant of the function f, see paper.
    %            If lambda is a vector, the result is a vector of the same
    %            shape and the values are independent element-wise results.
    %    xi_underline: The minimal value a verticy of S should have.
    %                  Shape (n + m x 1).
    %    xi_bar: The maximal value a verticy of S should have.
    %            Shape (n + m x 1).
    % 
    % Returns:
    %    sigma_upper_bound: The case (ii) upper bound of sigma.
    n_plus_m = size(xi_underline, 1);
    delta = norm(xi_bar - xi_underline);
    delta_s = sqrt(n_plus_m/(2*(n_plus_m + 1)))*delta;
    sigma_upper_bound = lambda * delta_s;
end
