function f = ode_euler(f_prime,x,f_o)
% x specifies the time points for which f should be evaluated
% f_o is the initial condition

% delta_x is the difference between two successive x values
delta_x = x(2)-x(1);

% determine how many points to approximate
l_x = length(x);

% initialize f
f = zeros(1,l_x);

% implement the equation of interest
for i = 1:(l_x-1)
    f(i+1) = feval(@f_prime,x(i),f_o);
end