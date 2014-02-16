% Parameters
% (- d_t + b * d_xx + a * d_x + c) u = 0

% Model Parameters
r = 0.05;
sigma = 1;

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
h = 1 / (m - 1); % 
domain_size = barrier;
grid_vector = linspace(0, domain_size, m);
num_el = numel(grid_vector);

%Time advancement parameters
initial_time = 0;
end_time = 1;

k = h; % k = O(h), makes sense for Crank-Nicolson
time_vector = initial_time:k:end_time;

%Boundary conditions
pattern_grid = zeros(size(grid_vector));
pattern_grid([1 end]) = 1;

pattern = logical(pattern_grid);

%We'll use the (discounted) initial conditions for the boundary. See Premia

%Initial conditions (Rainbow option)
initial_conditions_vector = max(grid_vector - K, 0)';

%Finite difference approximations
v = ones(num_el,1);
id_mat = speye(num_el);
tridiag = spdiags([v -2*v v], [-1 0 1], m, m); % D^2
centered_differences = spdiags([-v zeros(num_el, 1) v], [-1 0 1], m, m); % D

% u_t + H(Du, Iu) = G(u, D^2 U)
% Iu = \int{ u(x + z, t) - u(x, t) \Gamma_delta(z) dz}
% normpdf(X, mu, sigma)
% Trapezoidal rule: (1/2, 1, ..., 1, 1/2)
trap_rule_mat = ones(num_el, num_el);
trap_rule_mat(:, [1, end]) *= 0.5;
gamma_mat = normpdf(repmat(grid_vector, m, 1), 0, delta);
int_mat = trap_rule_mat .* gamma_mat;
