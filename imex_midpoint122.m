function [time_vector, solution] = imex_midpoint122(F, G, start, step, end_time)

  A = [0 0; ...
       0 0.5; ...
       0 1];

  A_hat = [0 0; ...
           0.5 0; ...
           0 1];

  [time_vector, solution] = rk_imex_solver(F, G, start, step, end_time, A, A_hat);
end

%!test
%! step = 0.001;
%! time_vector = 0:(step):1;
%!
%! F = @(t, x) -1;
%! G = [-1];
%! handle = @(t, x) G * x + F(t, x);
%!
%! [t, y] = ode45(handle, time_vector, [0.5]);
%!
%! [tt, yy] = imex_euler111(F, G, [0.5], step, 1);
%!
%! assert(yy, y', step);
