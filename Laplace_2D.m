% 2D Laplace Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

source setup_2D.m

%Crank-Nicolson [hint: this will take forever...]
%discretization_mat = b * tridiag + a * centered_differences + c * id_mat + lambda * int_mat;
%discretization_mat = b * laplacian + a * grad + c * id_mat;
discretization_mat = b * laplacian;
U_pruned = zeros(nnz(not(pattern)));

lhs_pruned = discretization_mat(not(pattern), not(pattern));
rhs_pruned = discretization_mat(not(pattern), pattern) * initial_conditions_vector(pattern);

%constant_term = ((discretization_mat - eye(m^2, m^2) / k) * pattern_matrix * initial_conditions_vector)(logical(1 - pattern_vector));
constant_term = discretization_mat(not(pattern), pattern) * initial_conditions_vector(pattern);

U_pruned = lhs_pruned \ rhs_pruned;

axis([0 barrier 0 barrier 0 300])
mesh(grid_x(2:m-1, 2:m-1), grid_y(2:m-1, 2:m-1), reshape(U_pruned(:), m-2, m-2));
%print(sprintf("frame%i.png", i), "-dpng")
pause(100)

input("Press any key to continue ..."); 
