function u_h = fin_elt_1d(h)
% fin_elt_1d solves 1-d Poisson's equation with Dirichlet
% boundary conditions for the specified homework function
%
% uh = fin_elt_1d(h) uses a mesh with where h(i) = x_{i+1} - x_i

    N = length(h);
    F = zeros(N-1,1);
    A = zeros(N-1, N-1);
    grid = cumsum(h);
    
    for i = 1:N-1
        A(i,i) = 1/h(i)+1/h(i+1);
        if i < N-1
            A(i,i+1) = -1/h(i+1);
            A(i+1,i) = -1/h(i);
        end
    end
    
    for i = 1:N-1
        
        %composite gaussian quadrature
        integral = 0;
        for j = 1:N
            if j == 1
                a = 0;
                b = grid(1);
            elseif j == N
                a = grid(N-1);
                b = 1;
            else
                a = grid(j-1);
                b = grid(j);
            end
            integral = integral + ((b-a)/2)*second_deriv( (-1/sqrt(3))*(b-a)/2 + (a+b)/2 )*phi(h, i, (-1/sqrt(3))*(b-a)/2 + (a+b)/2 );
            integral = integral + ((b-a)/2)*second_deriv( (1/sqrt(3))*(b-a)/2 + (a+b)/2 )*phi(h, i, (1/sqrt(3))*(b-a)/2 + (a+b)/2 ); 
        end
        F(i)= integral;  
    end
    
    u_h = A\F;
    