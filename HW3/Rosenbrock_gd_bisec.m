function [opt_path, conv_path] = Rosenbrock_gd_bisec(x0, y0, alpha_UB, epsilon, Nmax)
    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        
        % call rosenbrock.m at (xn,yn) to get the current fn and gn
        [fn, gn] = rosenbrock(xn, yn);
        % store fn in nth entry of conv_path, and store xn,yn in nth column
        % of opt path
        opt_path(:,n) = [xn; yn];
        conv_path(:, n) = fn;
        % check to see if norm of gradient is below epsilon tolerance - if
        % so, break the for loop.
        if norm(gn) < epsilon
            break;
        end
        % normalize the gradient to get \hat{g}, the unit vector in the
        % gradient direction
        g_hat = gn / norm(gn, 2);
    
        %% use bisection algorithm to find step size
        % define h'(alpha) for current iteration as described in main file
        
        x_new = @(alpha) xn - alpha * g_hat(1);
        y_new = @(alpha) yn - alpha*g_hat(2);
        grad = @(x, y) [-2*(1-x)-4*x*(y - x.^2); 2*(y - (x.^2))];

        h_alpha = @(alpha) dot(g_hat, grad(x_new(alpha), y_new(alpha)));
        [arr, ~] = bisection(h_alpha, 0, alpha_UB, 10e-10, 40);
        alpha_opt = arr(end);
        % take gradient descent step with optimal step size alpha. Note
        % that you must use the normalized gradient here because that's
        % what we used to calculate the optimal step size
        xn = xn - alpha_opt*g_hat(1);
        yn = yn - alpha_opt*g_hat(2);
        
    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn, yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end
