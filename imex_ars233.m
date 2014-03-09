function [time_vector, solution] = imex_ars233(F, G, start, step, end_time)

  gamma = (3 + sqrt(3)) / 6;
  A = [0 0 0; ...
       0 gamma 0; ...
       0 (1 - 2 * gamma) gamma;
       0 0.5 0.5];
       
  A_hat = [0 0 0; ...
           gamma 0 0; ...
           (gamma - 1) (2 * gamma - 2) 0;
           0 0.5 0.5];

  [time_vector, solution] = rk_imex_solver(F, G, start, step, end_time, A, A_hat);
end

%!test
%! step = 0.01;
%! atol = step * 9;
%! time_vector = 0:(step):1;
%!
%! F = @(t, x) [-1; -1];
%! G = [-1 0; 0 -1];
%! handle = @(t, x) G * x + F(t, x);
%!
%! [t, y] = ode45(handle, time_vector, [0.5]);
%!
%! [tt, yy] = imex_ars233(F, G, [0.5; 0.5], step, 1);
%!
%! assert(yy, y', atol);
