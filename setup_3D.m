% Parameters
% (- d_t + b * d_xx + a * d_x + c) u = 0

% Model Parameters
r = 0.05; %r = 0.05;
sigma = 0.3;

lambda = 1;
delta = 0.2;
E = exp(delta^2); % Moment generating function?

K = 20;
S_max = 60;
S_min = 1;

% PDE parameters, after a change of value
a = - (r - lambda * delta - sigma^2 /2);
b = sigma^2 / 2;
c = r - lambda;

if (abs(a) > abs(b))
  disp("Convezione dominante")
end

% Grid parameters
m = 61; % Number of grid points
h = (log(S_max) - log(S_min)) / (m - 1);

% x = log(S)
grid_vector = linspace(log(S_min), log(S_max), m);
[grid_x, grid_y, grid_z] = meshgrid(grid_vector);
num_el = numel(grid_x);

%Time advancement parameters
initial_time = 0;
end_time = 0.25;

k = h * 0.5; % k = O(h)
time_vector = initial_time:k:end_time;

%Boundary conditions
pattern_grid = zeros(size(grid_x));
pattern_grid([1 end], :, :) = 1;
pattern_grid(:, [1 end], :) = 1;
pattern_grid(:, :, [1 end]) = 1;

pattern = logical(reshape(pattern_grid, num_el, 1));

%We'll use the (discounted) initial conditions for the boundary. See Premia

%Initial conditions (Rainbow option)
initial_conditions_matrix = max(max(max(exp(grid_x), exp(grid_y)) - K, exp(grid_z)), 0);
initial_conditions_vector = reshape(initial_conditions_matrix, num_el, 1);

%Finite difference approximations
v = ones(num_el,1);
id_mat = speye(num_el);
tridiag = spdiags([v -2*v v], [-1 0 1], m, m); % D^2
centered_differences = spdiags([-v zeros(num_el, 1) v], [-1 0 1], m, m); % D

laplacian = sparse(num_el, num_el);
laplacian = kron(kron(tridiag, speye(m, m)), speye(m, m)) ...
          + kron(kron(speye(m, m), tridiag), speye(m, m)) ...
          + kron(kron(speye(m, m), speye(m, m)), tridiag);

mixed_derivatives = sparse(num_el, num_el);
mixed_derivatives = kron(kron(centered_differences, centered_differences), speye(m, m)) ...
                  + kron(kron(speye(m, m), centered_differences), centered_differences) ...
                  + kron(kron(centered_differences, speye(m, m)), centered_differences);
     
grad = sparse(num_el, num_el);
grad = kron(kron(centered_differences, speye(m, m)), speye(m, m)) ...
     + kron(kron(speye(m, m), centered_differences), speye(m, m)) ...
     + kron(kron(speye(m, m), speye(m, m)), centered_differences);

% u_t + H(Du, Iu) = G(u, D^2 U)
% Iu = \int{ u(x + z, t) - u(x, t) \Gamma_delta(z) dz}
% normpdf(X, mu, sigma)
% Trapezoidal rule: (1/2, 1, ..., 1, 1/2)
%trap_rule_mat = ones(1, num_el);
%trap_rule_mat(1, pattern) *= 0.5;
%FIXME the integral term is wrong. Build with the tensor product.
%gamma_vec = normpdf(reshape(grid_x, 1, num_el), 0, delta);
%gamma_vec = gamma_vec .* normpdf(reshape(grid_y, num_el, 1), 0, delta)';
%gamma_vec = gamma_vec .* normpdf(reshape(grid_z, num_el, 1), 0, delta)';
%gamma_vec(1, pattern) *= 0.5;

gamma_mat = ones(m, m);
gamma_mat(:, [1 end]) = gamma_mat(:, [1 end]) * 0.5; % Matlab dislikes *=

% We apply the change variables y = x + z, the convolution
%   \int u(x + z) Gamma(z) dz = \int u(y) Gamma(y - x) dy
% This means that the integral term in x is equal to the convolution of the
% solution with a gaussian random variable centered in x
gamma_mat = gamma_mat .* normpdf(grid_x(:, :, 1), grid_x(:, :, 1)', delta);
gamma_mat = sparse(gamma_mat);

int_mat = kron(kron(gamma_mat, speye(m, m)), speye(m, m));
int_mat = int_mat + kron(kron(speye(m, m), gamma_mat), speye(m, m));
int_mat = int_mat + kron(kron(speye(m, m), speye(m, m)), gamma_mat);

disp("Setup finished")
