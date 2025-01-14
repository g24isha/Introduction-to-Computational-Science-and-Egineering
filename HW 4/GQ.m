function I = GQ(func,d,n)
%GQ Compute the Gaussian quadrature rule for region Omega=[0,1]^d
%
%CALL:  I = GQ(func,d,n)
%  I = Gaussian quadrature approximation
%  func = function of interest
%  d = dimension of the region
%  n = number of abscissas along one dimension
assert(d==1|d==2,"Invalid dimension")

[x,w] = qrule(n);
x = 1/2*(x+1);
w = 1/2*w;

if d==1
    I = 0;
    for i = 1:length(w)
        I = I + w(i)*func(x(i));
    end
else
    I = 0;
    % [xx,yy] = meshgrid(x,x);
    % ww = w.*w';
    for i = 1:length(w)
        for j = 1:length(w)
            I = I + w(i)*w(j)*func(x(i), x(j));
        end
    end
end

end