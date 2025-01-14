function Iquad = quadrature_threepts(f,a,b,N)
Iquad = 0;
x = linspace(a,b,N);
h = x(2)-x(1);
fval = f(x);
for i = 1:(length(x)-1)/2
    Iquad = Iquad + threept_rule(h,fval(2*i-1),fval(2*i),fval(2*i+1)) ;
end
end


function I = threept_rule(h,f0,f1,f2)

I = 2*h*((5/12)*f0 + (2/3)*f1 + (-1/12)*f2);

end

