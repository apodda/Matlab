% 2D Parabolic Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

source setup.m

%Crank-Nicolson [hint: this will take forever...]
U_pruned = zeros(nnz(not(pattern)), length(time_vector));

discretization = b / (h^2) * laplacian + a / h * 0.5 * grad + c * speye(num_el, num_el);

lhs = (id_mat - 0.5 * k * discretization);
rhs = (id_mat + 0.5 * k * discretization);

lhs_pruned = lhs(not(pattern), not(pattern));
rhs_pruned = rhs(not(pattern), not(pattern));

constant_term = k * discretization(not(pattern), pattern) * initial_conditions_vector(pattern);

U_pruned(:,1) = initial_conditions_vector(not(pattern));

for i=2:length(time_vector)
  U_pruned(:, i) = bicgstab(lhs_pruned, (rhs_pruned * U_pruned(:, i-1)) + constant_term);
end

U = zeros(length(pattern), length(time_vector));

for i=1:length(time_vector)
  U(:, i) = initial_conditions_vector;
end

U(not(pattern), :) = U_pruned;

for i=1:length(time_vector)
  mesh(grid_x, grid_y, reshape(U(:, i), m, m));
  axis([-1 barrier -1 barrier 0 250])
  disp(sprintf("Graph at time %f", time_vector(i)))
  %print(sprintf("frame%i.png", i), "-dpng")
  pause(time_vector(i))
end

input("Press any key to continue ..."); 
