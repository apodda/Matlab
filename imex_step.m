function result = imex_step(start, step, A, b, c, A_hat, b_hat, c_hat, F, G)
  % TODO Error checking
  % F should be an handle, G a matrix.
  % TODO allow all combinations
  
  result = zeros(size(start));
  
  s = numel(c);
  tmp = zeros(numel(start), 1);
  
  explicit_stages = zeros(numel(start), s+1);
  implicit_stages = zeros(numel(start), s);
  
  explicit_stages(:, 1) = F(start);
  
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
    implicit_stages(:, i) = bicgstab(lhs, G * tmp, 10^-4, 60, [], [], start);
    
    % Solve X_expl = F * (tmp + a_ii * X_impl) + const
    explicit_stages(:, i+1) = F(tmp + step * A(i, i) * implicit_stages(:, i));
  end
  
  result = start + step * (implicit_stages * b' + explicit_stages * b_hat');
end

%!xtest
%! step = 0.1;
%! time_vector = 0:(step):1;
%!
%! A = [1];
%! b = [1];
%! c = [1];
%!
%! A_hat = [0 0; 1 0];
%! b_hat = [1 0];
%! c_hat = [0 1];
%!
%! handle = @(t, x) -x - 1;
%! F = @(x) -1;
%! G = [-1];
%!
%! [t, y] = ode45(handle, time_vector, [0.5]);
%!
%! v = y';
%! for i=3:length(v)
%!   v(i) = imex_step(v(i-1), step, A, b, c, A_hat, b_hat, c_hat, F, G);
%! end
%!
%! assert(v, y', step); 
