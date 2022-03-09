function [sigma_upper_bound] = compute_sigma_upper_bound(lambda, ...
                                                         xi_underline, ...
                                                         xi_bar)
    % Paper "Mesh-Based Affine Abstraction of Nonlinear Systems with 
    % Tighter Bounds", Proposition 1, case (ii).
    % 
    % Args:
    %    lambda: Lipschitz constant of the function f, see paper.
    %    xi_underline: The minimal value a verticy of S should have.
    %                  Shape (n + m x 1).
    %    xi_bar: The maximal value a verticy of S should have.
    %            Shape (n + m x 1).
    % 
    % Returns:
    %    sigma_upper_bound: The case (ii) upper bound of sigma.
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
