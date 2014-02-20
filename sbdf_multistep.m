function result = sbdf_multistep(start, step, F, G)
  % F handle (explicit), G matrix (implicit)
  %TODO Error checking
  result = zeros(size(start));
  
  % Solve 3 * U_n+1 - 4 * U_n + U_n-1 = 2k * (2 * F * U_n - F * U_n-1 + G * U^n+1)
  % => (3 * Id - 2k * G) U_n+1 = 4 * (Id + k * F) * U_n - (Id + 2 * k * F) U_n-1
  lhs = @(x) 3 * x - 2 * step * G * x;
  rhs = 4 * (start(:, 2) + step * F(start(:, 2))) - start(:, 1) - 2 * step * F(start(:, 1));
  %FIXME Do not hardcode tolerances
  result = bicgstab(lhs, rhs, 10^-4, 60, [], [], start(:, 2));
end

%!xtest
%! step = 0.1;
%! time_vector = 0:(step):1;
%!
%! handle = @(t, x) -x - 1;
%! F = @(x) -1;
%! G = [-1];
%!
%! [t, y] = ode45(handle, time_vector, [0.5]);
%!
%! v = y';
%! for i=3:length(v)
%!   v(i) = sbdf_multistep(v(i-2:i-1), 0.1, F, G);
%! end
%!
%! assert(v, y', step^2); 
