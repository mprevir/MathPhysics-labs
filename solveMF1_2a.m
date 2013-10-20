function solveMF1_2a( q_x, k, beta, gamma1, gamma2, f_x, l )
%solve following boundary value problem:
%  k*u''(x)-q(x)*u(x) = -f(x), 0<x<l
%  -k*u'(0)+beta*u(0)=gamma1, u(l)=gamma2
%
%usage:
%   solveMF1_2a( q_x, k, beta, gamma1, gamma2, f_x, l )
    
    h = l/500;
    x = 0:h:l;

    N = length(x);

    D = zeros(1, N); %main diagonal of matrix
    D(1) = k/h + beta;
    for i = 2:N-1
        D(i) = -2*k/(h^2) - q_x(x(i));
    end;
    D(N) = 1;

    B = zeros(1, N); %diagonal below the main
    for i = 2:N-1
        B(i) = k/(h^2);
    end;
    B(N) = 0;

    A = zeros(1, N); %diagonal above the main
    A(1) = -k/h;
    for i = 2:N-1
        A(i) = k/(h^2);
    end;
    
    %Tridiagonal matrix algorithm (Metod Progonki)
    % vector of free members: (gamma1, f_x(x(2)), f_x(x(3)),...,f_x(x(N-1)), gamma2)
    f = zeros(1,N);
    f(1) = gamma1;
    for i = 2:N-1
        f(i) = -f_x(x(i));
    end;
    f(N) = gamma2;
   
    a = zeros(1,N);
    b = zeros(1,N);
    a(1) = -A(1)/D(1);
    b(1) = f(1)/D(1); 
    for i = 2:N-1
        a(i+1) = -A(i) / (D(i) + B(i)*a(i));
        b(i+1) = (-B(i)*b(i)+f(i)) / (D(i)+B(i)*a(i));
    end;
    
    y = zeros(1,N);
    y(N) = (-B(N)*b(N)+f(N)) / (D(N)+B(N)*a(N));

    for i=N-1:-1:1
        y(i) = a(i+1)*y(i+1)+b(i+1);
    end;    

    plot( x, y );
end
