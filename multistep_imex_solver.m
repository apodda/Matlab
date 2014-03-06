function [time_vector, solution] = multistep_imex_solver(F, G, start, step, end_time, a, b, c)
  %TODO: error checking
  s = size(a, 2) - 1;
  assert(s <= size(start, 1), "Error: too many steps provided as a starting point");
  time_vector = 0:step:end_time;
  
  if(length(time_vector) < s)
    solution = start;
    error("Warning: Not enough time to inizialize the multistep method")
    return
  else
    solution = zeros(size(start, 1), length(time_vector));
    [rr, cc] = size(start);
    solution(1:rr, 1:cc) = start;
  end
  
  %Preallocate some space
  tmp = zeros(size(start, 1), 1);

  for jj=s+1:length(time_vector)
    tmp = step * F(time_vector(jj-s:jj-1), solution(:, jj-s:jj-1)) * b(1:s)' ...
        + step * (G * solution(:, jj-s:jj-1) * c(1:s)') ...
        - solution(:, jj-s:jj-1) * a(1:s)';
    if issparse(G)
      lhs = @(x) a(s+1) * x - step * c(s+1) * (G * x);
      solution(:, jj) = bicgstab(lhs, tmp, 10^-4, 60, [], [], solution(:, jj-1));
    else
      lhs = speye(size(G)) - step * c(s+1) * G;
      solution(:, jj) = lhs \ tmp;
    end
  end
end

%!test
%! a = [-1 1];
%! b = [1 0];
%! c = [0 1];
%!
%! step = 0.05;
%! time_vector = 0:(step):1;
%!
%! F = @(t, x) ones(size(x)) * t';
%! G = [-1 0; 0 -1];
%! handle = @(t, x) G * x + F(t, x);
%!
%! [t, y] = ode45(handle, time_vector, [0.5 0.5]);
%!
%! [tt, yy] = multistep_imex_solver(F, G, [0.5; 0.5], step, 1, a, b, c);
%!
%! assert(yy, y', step); 
