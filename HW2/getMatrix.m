function A = getMatrix(N, L)
    dx = L/N;
    vec1 = ones(1, N-2);
    vec2 = (-2)*ones(1, N-1);
    vec3 = (1)*ones(1, N-2);
    A2 = diag(vec1, 1) + diag(vec2) + diag(vec3, -1);
    %a is 0 and D is kappa, but that is ignored
    A = (A2/((dx)^2));
end