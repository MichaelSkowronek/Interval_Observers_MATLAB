% The illustrative example from "Interval Observers for Simultaneous 
% State and Model Estimation of Partially Known Nonlinear Systems" chapter
% V.
clear;
init_paths;

K = 250;
n = 2;
p = 1;
l = 3;
n_w = 3;
m = 1;
delta_t = 0.01;  % Forward Euler step size.
f = generate_f(@f_dot, delta_t);
h = generate_h(@h_dot, delta_t);
g = generate_g();
u_k = 0;
u = ones(K + 1, 1) * u_k;
v_bar = [0.1;
         0.1;
         0.1];
v_underline = - v_bar;
w_bar = v_bar;
w_underline = - v_bar;
x_0_bar = [0;
           0.6];
x_0_underline = [-0.35;
                 -0.1];
% TODO: What did they use here? [-10, 10] at least holds with respect to 
%       the example figures
d_0_underline = - 10;
d_0_bar = 10;
x_0 = [0;
       0.55];  % Some initial state within the bound similar to figures
d_0 = 0.05;  % Some initial state within the bound similar to figures
num_grid_points_per_dim_f_domain = 3;
num_grid_points_per_dim_g_domain = 3;
num_grid_points_per_dim_h_domain = 3;
num_interval_optimizations = 10;

% TODO: Calculate the next part myself
epsilon = 0.0001;
A_f = [0.994, -0.01, 1 - epsilon;
       0.009, 0.9965, - epsilon];  % Partial derivative lower bound
B_f = [1.006, -0.0065, 1 + epsilon;
       0.016, 1, epsilon];  % Partial derivative upper bound
C_f = [0, 0, 0;
       0, 0, 0];
f_d = generate_f_d(f, A_f, B_f, C_f, n, p, m, n_w);

%% Calculate Lipschitz constants
lambda_f_1 = calculate_lambda_f_1(delta_t, x_0_underline, x_0_bar);
lambda_f_2 = calculate_lambda_f_2(delta_t, x_0_underline, x_0_bar);
L_f = [lambda_f_1;
       lambda_f_2];
% grad_g_1 = [1; 0; 0; 0; 1; 0; 0], norm_grad_g_1 = sqrt(2)
% grad_g_2 = [0; 1; 0; 0; 0; 1; 0], norm_grad_g_2 = sqrt(2)
% grad_g_3 = [0; 0; cos(d); 0; 0; 0; 1], sq_norm_grad_g_3 = cos(d)^2 + 1
% We drop the d-interval constraint and get max(norm_grad_g_3) = sqrt(2).
L_g = [sqrt(2);
       sqrt(2);
       sqrt(2)];
% Gradient of h: grad_h = [- delta_t * 1/10 * sin(x_1);
%                          - delta_t * 1/10 * cos(x_2);
%                          1;
%                          0;
%                          0;
%                          0;
%                          delta_t]
% Squared norm of gradient of h:
%   sq_grad_h = delta_t^2 * 1/100 * (sin(x_1)^2 + cos(x_2)^2) +
%               delta_t^2 + 1
% We drop the x-interval constraint and find the maximum at 
% x = [Pi/2 + k_1 * Pi;
%      k_2 * Pi] with k_1, k_2 some integer.
L_h = sqrt(delta_t^2 * 1/100 * 2 + delta_t^2 + 1);

%% Simulation to get y
delta_t_sim = 0.00001;
sim_num_steps_factor = 1000;  % delta_t / delta_t_sim
f_sim = generate_f(@f_dot, delta_t_sim);
h_sim = generate_h(@h_dot, delta_t_sim);
num_simulation_steps = K * sim_num_steps_factor;
u_sim = ones(num_simulation_steps  + 1, 1) * u_k;
[~, ~, y_sim] = simulate_system(f_sim, h_sim, g, n, p, l, n_w, x_0, ...
                                d_0, u_sim, w_underline, w_bar, ...
                                v_underline, v_bar, num_simulation_steps);
y = y_sim(1:sim_num_steps_factor:end, :);

%%
smio(f, f_d, g, y, u, x_0_underline, ...
     x_0_bar, d_0_underline, d_0_bar, w_underline, w_bar, ...
     v_underline, v_bar, L_f, L_g, L_h, ...
     num_grid_points_per_dim_f_domain, ...
     num_grid_points_per_dim_g_domain, ...
     num_grid_points_per_dim_h_domain, K, num_interval_optimizations, ...
     'plot_global_affine_abstraction_f', false, 'x_spacing', 0.001, ...
     'y_spacing', 0.001);


function f_handle = generate_f(f_dot, delta_t)
    function x_k_plus_one = f(xi_k)
        x_k = xi_k(1:2);
        x_dot_k = f_dot(xi_k);
        x_k_plus_one = x_k + delta_t * x_dot_k;
    end
    f_handle = @f;
end


function x_dot = f_dot(xi)
x = xi(1:2);
d = xi(3);
u = xi(4);
w = xi(5:7);
x_dot = zeros(2, 1);
x_dot(1) = - x(1) * x(2) - x(2) + u + d + w(1);
x_dot(2) = x(1) * x(2) + x(1) + w(2);
end


function f_1_handle = generate_f_1(f_1_dot, delta_t)
    function x_1_k_plus_one = f_1(xi_k)
        x_1_k = xi_k(1);
        x_1_dot_k = f_1_dot(xi_k);
        x_1_k_plus_one = x_1_k + delta_t * x_1_dot_k;
    end
    f_1_handle = @f_1;
end


function x_1_dot = f_1_dot(xi)
x = xi(1:2);
d = xi(3);
u = xi(4);
w = xi(5:7);
x_1_dot = - x(1) * x(2) - x(2) + u + d + w(1);
end


function f_2_handle = generate_f_2(f_2_dot, delta_t)
    function x_2_k_plus_one = f_2(xi_k)
        x_2_k = xi_k(2);
        x_2_dot_k = f_2_dot(xi_k);
        x_2_k_plus_one = x_2_k + delta_t * x_2_dot_k;
    end
    f_2_handle = @f_2;
end


function x_2_dot = f_2_dot(xi)
x = xi(1:2);
w = xi(5:7);
x_2_dot = x(1) * x(2) + x(1) + w(2);
end


function f_d_handle = generate_f_d(f, A_f, B_f, C_f, n, p, m, n_w)
    function sol = f_d(xi_x, xi_y)
        % (1) from "Interval Observers for Simultaneous State and Model 
        % Estimation of Partially Known Nonlinear Systems"
        %
        % Args:
        %   xi_x: Vector of size (n + p + m + n_w) which consists of x, d, 
        %         u, w.
        %   xi_y: Vector of size (n + p + m + n_w) which consists of x, d,
        %         u, w.
        % TODO:
        %   Why are u and w not considered? We need to choose u and w to
        %   evaluate f. As u is always the same for this example, just
        %   choosing the u might be ok. w = 0 also might be ok.
        x = xi_x(1:n + p);
        y = xi_y(1:n + p);
        u = xi_x(n + p + 1:n + p + m);
        w = zeros(n_w, 1);
        m_d = n;  % Note: m_d is named m within literature
        n_d = n + p;  % Note: n_d is named n within literature
        sol = zeros(m_d, 1);
        for i = 1:m_d
            z = zeros(n_d, 1);
            for j = 1:n_d
                a_i_j = A_f(i, j);
                b_i_j = B_f(i, j);
                if a_i_j >= 0 || ...
                   (a_i_j <= 0 && b_i_j >= 0 && abs(a_i_j) <= abs(b_i_j))
                    % Case 1 or 2 from "On sufficient conditions for mixed 
                    % monotonicity", (11).
                    z(j) = x(j);
                else
                    % Case 3 or 4 from "On sufficient conditions for mixed 
                    % monotonicity", (11).
                    z(j) = y(j);
                end
            end
            f_value = f([z;
                         u;
                         w]);  % TODO: Split f() into f_1(), f_2(), ...
            sol(i) = f_value(i) + C_f(i, :) * (x - y); 
        end
    end
    f_d_handle = @f_d;
end


function h_handle = generate_h(h_dot, delta_t)
    function d_k_plus_one = h(xi_k)
        d_k = xi_k(3);
        d_k_plus_one = d_k + delta_t * h_dot(xi_k);
    end
    h_handle = @h;
end


function d_dot = h_dot(xi)
x = xi(1:2);
w = xi(5:7);
d_dot = 0.1 * (cos(x(1)) - sin(x(2))) + w(3);
end


function g_handle = generate_g()
    function y = g(xi_g)
        % g([x; d; u; v]) with xi_g = [x; d; u; v].
        x = xi_g(1:2);
        d = xi_g(3);
        v = xi_g(5:7);
        y = zeros(3, 1);
        y(1) = x(1) + v(1);
        y(2) = x(2) + v(2);
        y(3) = sin(d) + v(3);
    end
    g_handle = @g;
end


function [x, d, y] = simulate_system(f, h, g, n, p, l, n_w, x_0, d_0, ...
                                     u, w_underline, w_bar, ...
                                     v_underline, v_bar, ...
                                     num_simulation_steps)
    % Simulates the system (f, h, g) for num_simulation_steps.
    %
    % Args:
    %   u: A matrix of size (num_simulation_steps + 1 x m) of input 
    %      vectors. The first entry is the input of step 0.
    %   K: The number of simulation steps.
    % Returns:
    %   x: Size (num_simulation steps + 1 x n) starting with x_0.
    %   d: Size (num_simulation steps + 1 x p) starting with d_0.
    %   y: Size (num_simulation steps + 1 x l) starting with y_0.
    % TODO:
    %   Consider other distributions
    w_samples = rand(num_simulation_steps + 1, n_w);
    w = w_underline' + w_samples * (w_bar - w_underline);
    v = zeros(num_simulation_steps + 1, l);
    %v_samples = rand(simulation_steps + 1, l);
    %v = v_underline' + v_samples * (v_bar - v_underline);
    x = zeros(num_simulation_steps + 1, n);
    d = zeros(num_simulation_steps + 1, p);
    y = zeros(num_simulation_steps + 1, l);
    x(1, :) = x_0';
    d(1, :) = d_0';
    
    u_0 = u(1, :)';
    v_0 = v(1, :)';
    xi_0_g = [x_0;
              d_0;
              u_0;
              v_0];
    y_0 = g(xi_0_g);
    y(1, :) = y_0';
    for k = 1:num_simulation_steps
        x_k_minus_one = x(k, :)';
        d_k_minus_one = d(k, :)';
        u_k_minus_one = u(k, :)';
        w_k_minus_one = w(k, :)';
        xi_k_minus_one = [x_k_minus_one;
                          d_k_minus_one;
                          u_k_minus_one;
                          w_k_minus_one];
        x_k = f(xi_k_minus_one);
        d_k = h(xi_k_minus_one);
        x(k + 1, :) = x_k';
        d(k + 1, :) = d_k';
        
        u_k = u(k + 1, :)';
        v_k = v(k + 1, :)';
        xi_k_g = [x_k;
                  d_k;
                  u_k;
                  v_k];
        y_k = g(xi_k_g);
        y(k + 1, :) = y_k';
    end
end


function lambda_f_1 = calculate_lambda_f_1(delta_t, x_0_underline, x_0_bar)
    % Calculate Lipschitz constant of f_1.
    %
    % Notation follows the MATLAB quadprog notation.
    % Gradient of f_1: grad_f_1 = [1 - delta_t * x_2;
    %                              - delta_t * (x_1 + 1);
    %                              delta_t;
    %                              delta_t;
    %                              delta_t;
    %                              0;
    %                              0]
    % Squared norm of gradient of f_1: 
    %   sq_grad_f_1 = delta_t^2 * x_1^2 + 2 * delta_t^2 * x_1 + 
    %                 delta_t^2 * x_2^2 - 2 * delta_t * x_2 +
    %                 4 * delta_t^2 + 1
    %
    % Note:
    %   - We need to negate H and f in order to maximize.
    %   - We only use x_1 and x_2 because the objective is independent of
    %     the remaining variables.
    H = - [2 * delta_t^2, 0;
           0            , 2 * delta_t^2];
    f = - [2 * delta_t^2;
           - 2 * delta_t];
    lb = x_0_underline;
    ub = x_0_bar;
    x = quadprog(H, f, [], [], [], [], lb, ub);
    x_1 = x(1);
    x_2 = x(2);
    lambda_f_1 = sqrt(delta_t^2 * x_1^2 + 2 * delta_t^2 * x_1 + ...
                      delta_t^2 * x_2^2 - 2 * delta_t * x_2 + ...
                      4 * delta_t^2 + 1);
end


function lambda_f_2 = calculate_lambda_f_2(delta_t, x_0_underline, x_0_bar)
    % Calculate Lipschitz constant of f_2.
    %
    % Notation follows the MATLAB quadprog notation.
    % Gradient of f_2: grad_f_2 = [delta_t * (x_2 + 1);
    %                              1 + delta_t * x_1;
    %                              0;
    %                              0;
    %                              0;
    %                              delta_t;
    %                              0]
    % Squared norm of gradient of f_2: 
    %   sq_grad_f_2 = delta_t^2 * x_1^2 + 2 * delta_t * x_1 +  
    %                 delta_t^2 * x_2^2 + 2 * delta_t^2 * x_2 + 
    %                 2 * delta_t^2
    %
    % Note:
    %   - We need to negate H and f in order to maximize.
    %   - We only use x_1 and x_2 because the objective is independent of
    %     the remaining variables.
    H = - [2 * delta_t^2, 0;
           0            , 2 * delta_t^2];
    f = - [2 * delta_t;
           2 * delta_t^2];
    lb = x_0_underline;
    ub = x_0_bar;
    x = quadprog(H, f, [], [], [], [], lb, ub);
    x_1 = x(1);
    x_2 = x(2);
    lambda_f_2 = sqrt(delta_t^2 * x_1^2 + 2 * delta_t * x_1 + ...
                      delta_t^2 * x_2^2 + 2 * delta_t^2 * x_2 + ...
                      2 * delta_t^2);
end
