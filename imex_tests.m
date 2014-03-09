%!shared tol, step, start, end_time, F, G, y
%! more off
%! tol = 0.1;
%!
%! step = 0.05;
%! start = [0.5 0.5];
%! end_time = 1;
%! time_vector = 0:(step):end_time;
%!
%! %F = @(t, x) ones(size(x)) .* t;
%! F = @(t, x) zeros(size(x));
%! G = [-1 0; 0 -1];
%! handle = @(t, x) G * x + F(t, x);
%!
%! [t, y] = ode45(handle, time_vector, start);


%% Imex Multistep tests

%!test
%! [tt, yy] = imex_cnab(F, G, start, step, end_time);
%! assert(yy, y', tol); 

%!test
%! [tt, yy] = imex_mcnab(F, G, start, step, end_time);
%! assert(yy, y', tol); 

%!xtest
%! [tt, yy] = imex_cnlf(F, G, start, step, end_time);
%! assert(yy, y', tol); 

%!test
%! [tt, yy] = imex_sbdf2(F, G, start, step, end_time);
%! assert(yy, y', tol); 

%% Imex Runge-Kutta tests

%!test
%! [tt, yy] = imex_euler111(F, G, start, step, end_time);
%! assert(yy, y', tol);

%!test
%! [tt, yy] = imex_midpoint122(F, G, start, step, end_time);
%! assert(yy, y', tol);

%!test
%! [tt, yy] = imex_ars233(F, G, start, step, end_time);
%! assert(yy, y', tol);
