function result = sbdf_multistep(start, step, F, G, f)
  % F explicit, G implicit, f known term
  %TODO Error checking
  result = zeros(size(start));
  
  % Solve 3 * U_n+1 - 4 * U_n + U_n-1 = 2k * (2 * F * U_n - F * U_n-1 + G * U^n+1)
  % => (3 * Id - 2k * G) U_n+1 = 4 * (Id + k * F) * U_n - (Id + 2 * k * F) U_n-1
  lhs = @(x) 3 * x - 2 * step * G * x;
  rhs = 4 * (start(:, 2) + step * F * start(:, 2) + step * f) - ...
        start(:, 1) - step * F * start(:, 1) - 2 * step * f;
  result = bicgstab(lhs, rhs, 10^-6, 60, [], [], start(:, 2));
end
