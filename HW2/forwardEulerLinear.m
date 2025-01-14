function uT = forwardEulerLinear(u0,kappa,A,dt,T)
    % implement forward Euler here and return just the solution at the
    % final time.
    N = length(u0);
    K = ceil(T/dt);
    u = zeros(N, K);
    u(:, 1) = u0 + kappa*dt*A*u0;
    for i=2:K
        u(:, i) = u(:, i -1) + kappa*dt*A*u(:, i-1);
    end
    uT = u(:, end);
end
