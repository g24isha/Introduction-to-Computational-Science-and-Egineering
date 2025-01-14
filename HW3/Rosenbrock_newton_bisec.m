function [opt_path, conv_path] = Rosenbrock_newton_bisec(x0, y0, alpha_UB, epsilon, Nmax)

    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        % call rosenbrock.m at (xn,yn) to get the current fn, gn, and Hn
        [fn, gn, Hn] = rosenbrock(xn, yn);
        
        % store fn in nth entry of conv_path, and store xn,yn in nth column
        % of opt path
        opt_path(:,n) = [xn; yn];
        conv_path(:, n) = fn;
        % check to see if norm of gradient is below epsilon tolerance - if
        % so, break the for loop.
        if norm(gn, 2) < epsilon
            break;
        end
        
        %% use bisection algorithm to find step size
        % define \hat{g}, the Newton update direction, as defined in the
        % main file - note that this is NOT the same as \hat{g} for
        % gradient descent
        d = Hn\gn;
        g_hat = d/norm(d, 2);
        
        % check if \hat{g} is a descent direction. If not, change the
        % sign of \hat{g}
        gn_transpose = gn';
        if gn_transpose*g_hat < 0
            g_hat = (-1) * g_hat;
        end

        % define q(alpha) for current iteration as described in main file
        x_new = @(alpha) xn - alpha * g_hat(1);
        y_new = @(alpha) yn - alpha*g_hat(2);
        grad = @(x, y) [-2*(1-x)-4*x*(y - x.^2); 2*(y - (x.^2))];

        q_alpha = @(alpha) dot(g_hat, grad(x_new(alpha), y_new(alpha)));
        [arr, ~] = bisection(q_alpha, 0, alpha_UB, 10e-10, 40);
        alpha_opt = arr(end);
        % call bisection routine on q(alpha) to find optimal step size
        
        % take Newton step with optimal step size alpha (you must use the
        % normalized Newton direction here)
        xn = xn - alpha_opt*g_hat(1);
        yn = yn - alpha_opt*g_hat(2);

    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn,yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end