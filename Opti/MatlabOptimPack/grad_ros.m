function g = grad_ros(x);
% Compute the gradient of the rosenbrock function

g = zeros_(size(x));
g(1) = -2.0*(1.0-x(1)) - 400.0*(x(2)-x(1).^2)*x(1);
g(2) = 200.0*(x(2)-x(1).^2);
 