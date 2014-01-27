% 2D Laplace Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

% Model Parameters
r = 0.05;
sigma = 0.2;

lambda = 1;
delta = 1;

K = 100;
barrier = 300; %FIXME remove/rename

% PDE parameters, after a change of value
a = - (r - sigma^2 /2); %FIXME
b = sigma^2 / 2;
c = r;

% Grid parameters
m = 100; % Number of grid points
num_el = m^2;
h = 1 / (m - 1); % 
domain_size = barrier;
grid_vector = linspace(-domain_size, domain_size, m);
[grid_x, grid_y] = meshgrid(grid_vector);

%Boundary conditions
pattern_grid = zeros(size(grid_x));
pattern_grid([1 end], :) = 1;
pattern_grid(:, [1 end]) = 1;

pattern = logical(reshape(pattern_grid, num_el, 1));

%We'll use the (discounted) initial conditions for the boundary. See Premia

%Initial conditions (Rainbow option)
initial_conditions_matrix = max(max(grid_x, grid_y) - K, 0);
initial_conditions_vector = reshape(initial_conditions_matrix, num_el, 1);

%Finite difference approximations
v = ones(num_el,1);
id_mat = speye(num_el);
tridiag = spdiags([v -2*v v], [-1 0 1], m, m) * (1 / h^2); % D^2
centered_differences = spdiags([-v zeros(num_el, 1) v], [-1 0 1], m, m) * (1 / (2 * h)); % D

laplacian = sparse(num_el, num_el);
laplacian = kron(tridiag, speye(m, m)) + kron(speye(m, m), tridiag);

grad = sparse(num_el, num_el);
grad = kron(centered_differences, speye(m, m)) + kron(speye(m, m), centered_differences);

% u_t + H(Du, Iu) = G(u, D^2 U)
% Iu = \int{ u(x + z, t) - u(x, t) \Gamma_delta(z) dz}
% normpdf(X, mu, sigma)
% Trapezoidal rule: (1/2, 1, ..., 1, 1/2)
%trap_rule_mat = ones(m, m);
%trap_rule_mat(:, [1, end]) *= 0.5;
%gamma_mat = normpdf(repmat(grid_matrix, m, 1), 0, delta);
%int_mat = trap_rule_mat .* gamma_mat;

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
