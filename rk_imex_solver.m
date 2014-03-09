function [time_vector, solution] = rk_imex_solver(F, G, start, step, end_time, tableau_implicit, tableau_explicit)
  
  time_vector = 0:step:end_time;
  solution = zeros(length(start), length(time_vector));
  
  solution(:, 1) = start;
  
  % Extract the parameters from the tableau
  A = tableau_implicit(2:(end-1), 2:end);
  A_hat = tableau_explicit(1:(end-1), 1:end);

  b = tableau_implicit(end, 2:end);
  b_hat = tableau_explicit(end, 1:end);
  
  c = sum(A');
  c_hat = sum(A_hat');
  
  tmp = zeros(numel(start), 1);
  s = numel(c);
  
  explicit_stages = zeros(numel(start), s+1);
  implicit_stages = zeros(numel(start), s);
  
  for jj=2:length(time_vector)
    explicit_stages(:, 1) = F(time_vector(jj) + c_hat(1) * step, solution(:, jj-1));
  
    for ii=1:s
      if ii==1
        tmp = solution(:, jj-1) + step * (explicit_stages(:, 1:ii) * A_hat(ii+1, 1:ii)');
      else
        tmp = solution(:, jj-1) + step * (explicit_stages(:, 1:ii) * A_hat(ii+1, 1:ii)' + ...
                              implicit_stages(:, 1:(ii-1)) * A(ii, 1:(ii-1))' );
      end
      
      if issparse(G)
        % Solve X = G * (tmp + step * a_ii * X) => (I - step * a_ii * G) X = G * tmp
        lhs = @(x) x - step * A(ii, ii) * G * x;
        %FIXME Do not hardcode tolerances
        %[L, U] = luinc(G, 10^-1);
        L = [];
        U = [];
        implicit_stages(:, ii) = bicgstab(lhs, G * tmp, 10^-4, 60, L, U, solution(:, jj-1));
      else
        lhs = speye(length(start)) - step * A(ii, ii) * G;
        implicit_stages(:, ii) = lhs \ (G * tmp);
      end
      
      % Solve X_expl = F * (tmp + a_ii * X_impl) + const
      explicit_stages(:, ii+1) = F(time_vector(ii) + c_hat(s) * step, tmp + step * A(ii, ii) * implicit_stages(:, ii));
    end
    
    solution(:, jj) = solution(:, jj-1) + step * (implicit_stages * b' + explicit_stages * b_hat');
  end
end
