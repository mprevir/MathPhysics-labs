function solveMF1_1( f_x, l )
%solves the equation u''_x = -f_x; u(0)=u(l)=0  using Fourier's method
%
%usage:
%   sovleMF1_1( f_x, l )

    h = l/1500;
    x = 0:h:l;

    N = length(x)-1;    
    
    lambda = (4/(h^2)) * (sin(pi*(1:N-1)*h/(2*l))).^2;

    mu_x = sqrt(2/l) * sin(pi*(1:N-1)' * x(2:N)/l);    

    Frr = arrayfun(f_x, x(2:N)) * mu_x * h;

    C = Frr./lambda;
    
    y = C * mu_x;

    plot(x(1:N-1), y);
end
