% 2D Parabolic Equation
% (- d_t + b * d_xx + a * d_x + c) u = 0

source setup_3D.m

%IMEX, Forward-Backward Euler
% u_n = u_{n-1} + k(f(u_{n-1}) + g(u_n))
%implicit_disc = b / (h^2) * laplacian + c * speye(num_el, num_el);
%explicit_disc = a / h * 0.5 * grad;

%lhs = id_mat - k * implicit_disc;
%rhs = id_mat + k * explicit_disc;

%constant_term = k * (implicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern) + ...
%                     explicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern));

%lhs_pruned = lhs(not(pattern), not(pattern));
%rhs_pruned = rhs(not(pattern), not(pattern));

% Setup
discretization = b / (h^2) * laplacian + a / h * 0.5 * grad + c * speye(num_el, num_el);

discretization_pruned = discretization(not(pattern), not(pattern));

f_cn = discretization(not(pattern), pattern) * initial_conditions_vector(pattern);

%IMEX
implicit_disc = b / (h^2) * laplacian + c * speye(num_el, num_el) + 0.25 * mixed_derivatives / (h^2);
explicit_disc = a / h * 0.5 * grad;
%gamma = (3 + sqrt(3))/6;

%A = [gamma 0; (1 - 2 * gamma) gamma];
%b = [0.5 0.5];
%c = [gamma (1 - gamma)];

%A_hat = [0 0 0; gamma 0 0; (gamma - 1) 2*(1 - gamma) 0];
%b_hat = [0 0.5 0.5];
%c_hat = [0 c];

A = [1];
b = [1];
c = [1];

A_hat = [0 0; 1 0];
b_hat = [1 0];
c_hat = [0 1];

U_pruned = zeros(nnz(not(pattern)), length(time_vector));
U_pruned(:,1) = initial_conditions_vector(not(pattern));

f = implicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern) + ...
    explicit_disc(not(pattern), pattern) * initial_conditions_vector(pattern);

tic
%Bootstrapping with Forward-Backward Euler. This has a smoothing effect on
% initial conditions
U_pruned(:, 2) = imex_step(U_pruned(:, 1), k, A, b, c, A_hat, b_hat, c_hat, ...
                       explicit_disc(not(pattern), not(pattern)), ...
                       implicit_disc(not(pattern), not(pattern)), ...
                       f);

for i=3:length(time_vector)
  % 2nd order SBDF?
  U_pruned(:, i) = sbdf_multistep(U_pruned(:,i-2:i-1), k, ...
                                  explicit_disc(not(pattern), not(pattern)), ...
                                  implicit_disc(not(pattern), not(pattern)), f);
end
toc

U = zeros(length(pattern), length(time_vector));

for i=1:length(time_vector)
  U(:, i) = initial_conditions_vector;
end

U(not(pattern), :) = U_pruned;

for i=1:length(time_vector)
  mesh(grid_x(:, :, 1), grid_y(:, :, 1), reshape(U(:, i), m, m, m)(:, :, 20));
  %axis([-1 barrier -1 barrier 0 250])
  disp(sprintf("Graph at time %f", time_vector(i)))
  %print(sprintf("frame%i.png", i), "-dpng")
  pause(time_vector(i))
end

input("Press any key to continue ..."); 
