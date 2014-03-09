function [time_vector, solution] = imex_cnab(F, G, start, step, end_time)

  a = [0 -1 1];
  b = [-1/2 3/2 0];
  c = [0 1/2 1/2];

  % FIXME handle better the case of 
  if (step > end_time)
    disp("Warning: Couldn't reach multistep phase. Not enough steps")
    [time_vector, solution] = imex_euler111(F, G, start, end_time, end_time);
    return;
  end

  [tt, bootstrap] = imex_euler111(F, G, start, step, step);

  [time_vector, solution] = multistep_imex_solver(F, G, bootstrap, step, end_time, a, b, c);
end

%!test
%!
%! step = 0.05;
%! time_vector = 0:(step):1;
%!
%! F = @(t, x) ones(size(x)) .* t;
%! G = [-1 0; 0 -1];
%! handle = @(t, x) G * x + F(t, x);
%!
%! [t, y] = ode45(handle, time_vector, [0.5 0.5]);
%!
%! [tt, yy] = imex_cnab(F, G, [0.5; 0.5], step, 1);
%!
%! assert(yy, y', step); 
