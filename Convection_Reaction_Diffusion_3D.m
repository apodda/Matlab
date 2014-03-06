% 2D Parabolic Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

source setup_3D.m

%IMEX
implicit_disc = b / (h^2) * laplacian + c * speye(num_el, num_el) + 0.25 * sigma * sigma * 0.5 * mixed_derivatives / (h^2);
explicit_disc = a / h * 0.5 * grad;
integral_term = lambda * h^3 * gamma_mat(not(pattern));

iv = initial_conditions_vector(not(pattern)');

f = implicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern) + ...
    explicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern);

F = @(t, x) explicit_disc(not(pattern), not(pattern)) * x + exp(-r * t) * f; % + integral_term * x * ones(size(x));
G = implicit_disc(not(pattern), not(pattern));

%gamma = (3 + sqrt(3)) / 6;
%k = sqrt(4 * gamma - 1) * h / (3 * gamma * a);

tic
[time_vector U_pruned] = imex_midpoint122(F, G, iv, k, end_time);
toc

U = zeros(length(pattern), length(time_vector));

for i=1:length(time_vector)
  U(:, i) = initial_conditions_vector;
end

U(not(pattern), :) = U_pruned;

gui_mode("3d")
colormap summer
shading interp

for i=1:length(time_vector)
  surf(exp(grid_x(:, :, 1)), exp(grid_y(:, :, 1)), reshape(U(:, i), m, m, m)(:, :, 20));
  axis([S_min S_max S_min S_max 0 250])
  disp(sprintf("Graph at time %f", time_vector(i)))
  %print(sprintf("frame%i.png", i), "-dpng")
  pause(time_vector(i))
end

input("Press any key to continue ..."); 
