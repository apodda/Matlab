function result = imex_step(start, step, A, b, c, A_hat, b_hat, c_hat, F, G, f)
  %TODO Error checking
  
  result = zeros(size(start));
  
  s = numel(c);
  tmp = zeros(numel(start), 1);
  
  explicit_stages = zeros(numel(start), s+1);
  implicit_stages = zeros(numel(start), s);
  
  explicit_stages(:, 1) = F * start + f;
  
  for i=1:s
    if i==1
      tmp = start + step * (explicit_stages(:, 1:i) * A_hat(i+1, 1:i)');
    else
      tmp = start + step * (explicit_stages(:, 1:i) * A_hat(i+1, 1:i)' + ...
                            implicit_stages(:, 1:(i-1)) * A(i, 1:(i-1))' );
    end
    
    % Solve X = G * (tmp + step * a_ii * X) => (I - step * a_ii * G) X = G * tmp
    lhs = @(x) x - step * A(i, i) * G * x;
    %FIXME Do not hardcode tolerances
    implicit_stages(:, i) = bicgstab(lhs, G * tmp, 10^-6, 60, [], [], start);
    
    % Solve X_expl = F * (tmp + a_ii * X_impl) + const
    explicit_stages(:, i+1) = F * (tmp + step * A(i, i) * implicit_stages(:, i)) + f;
  end
  
  result = start + step * (implicit_stages * b' + explicit_stages * b_hat');
end
