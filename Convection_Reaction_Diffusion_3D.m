% 2D Parabolic Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

source setup_3D.m

%IMEX
implicit_disc = b / (h^2) * laplacian + c * speye(num_el, num_el) + 0.25 * sigma * sigma * 0.5 * mixed_derivatives / (h^2);
explicit_disc = a / h * 0.5 * grad + lambda * h^3 * int_mat;

iv = initial_conditions_vector(not(pattern)');

f = implicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern) + ...
    explicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern);

F = @(t, x) explicit_disc(not(pattern), not(pattern)) * x + f * exp(-r * t);
G = implicit_disc(not(pattern), not(pattern));

clear implicit_disc explicit_disc f grad laplacian mixed_derivatives int_mat gamma_mat;

%gamma = (3 + sqrt(3)) / 6;
%k = sqrt(4 * gamma - 1) * h / (3 * gamma * a);

tic
%[time_vector U_pruned] = imex_midpoint122(F, G, iv, k, end_time);
[time_vector, U_pruned] = imex_sbdf2(F, G, iv, k, end_time);
toc

clear F G;

U = initial_conditions_vector * exp(-r * time_vector);
assert(size(U) == [length(pattern), length(time_vector)]);

U(not(pattern), :) = U_pruned;

clear U_pruned;

graphics_toolkit("fltk");
gui_mode("3d");
colormap summer
shading interp

slice = floor(m/2);

for i=1:length(time_vector)
  %surf(exp(grid_x(:, :, 1)), exp(grid_y(:, :, 1)), reshape(U(:, i), m, m, m)(:, :, 20), 'EdgeColor','none','LineStyle','none');
  surf(exp(grid_x(:, :, 1)), exp(grid_y(:, :, 1)), reshape(U(:, i), m, m, m)(:, :, slice));
  axis([S_min S_max S_min S_max 0 S_max])
  disp(sprintf("Graph at time %f", time_vector(i)))
  %print(sprintf("frame%i.png", i), "-dpng")
  if i>1
    tt = time_vector(i) - time_vector(i-1);
  else
    tt = 0;
  end
  pause(tt * 10)
end

input("Press any key to continue ..."); 
