function uT = backwardEulerLinear(u0,kappa,A,dt,T)
    N = length(u0);
    K = ceil(T/dt);

    % implement backward Euler here and output just the solution at the
    % final time.
    N = length(u0);
    K = ceil(T/dt);
    vec1 = ones(1, N);
    I = diag(vec1);
    u = zeros(N, K);
    u(:, 1) = (inv(I - kappa*dt*A))*u0';
    for i=2:K
        u(:, i) = (inv(I - kappa*dt*A))*u(:, i-1);
    end
    uT = u(:, end);
end
