% Plotting the Euler spiral
fx = @ (x) cos(x.^2); 
fy = @ (x) sin(x.^2); 

X = @(t) quadrature_threepts(fx,0,t,1050); 
Y = @(t) quadrature_threepts(fy,0,t,1050); 

t = -40:0.01:40; 
Xvec = zeros(length(t),1); 
Yvec = zeros(length(t),1); 
for i = 1:length(t)
    Xvec(i) = X(t(i)); 
    Yvec(i) = Y(t(i)); 
end

plot(Xvec,Yvec); 
grid on; 
