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
