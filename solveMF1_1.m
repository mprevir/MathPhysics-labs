function solveMF1_1( f_x, l )
%solves the equation u''_x = -f_x using Fourier's method

    h = l/150;
    x = 0:h:l;

    N = length(x) - 1;    
    
    lambda = zeros(1, N-1); %eigenvalues
    for k = 2:N-1
        lambda(k) = ( 4/(h^2)) * (sin(pi*k*h/(2*l)))^2;
    end;

    mu_x = zeros(N-1, N-1); %eigenfunctions
    for k = 2:N-1
        for i = 2:N-1
            mu_x(k, i) = sqrt(2/l) * sin(pi*k*x(i)/l);
        end;
    end;

    Frr = zeros(1, N-1); %Fourier's decomposition coefficients
    for k = 2:N-1;
        for j = 2:N-1
            Frr(k) += f_x(x(j)) * mu_x(k,j) * h;
        end;
    end;

    C = Frr./lambda;

    y = zeros(1, N+1); %solution
    for i = 2:N-1
        for k = 2:N-1
            y(i) += C(k) * mu_x(k, i);
        end; 
    end;

    plot(x, y);
end
