function [uh, uh_prime] = u_heval(u_h, grid, x)    
% u_heval creates the piecewise linear solution to Poisson's equation
%
% [uh, uh_prime] = fin_elt_1d(u_h, grid, x) gives the evaluation of the 
% piecewise linear interpolant over the grid as well as the derivative at 
% a point x. 
% The grid corresponds with the approximated solution u_h, calculated from
% fin_elt_1d.

    %append values at the boundaries
    grid = [0; grid; 1];
    u_h = [0; u_h; 0];
    lower = 0; upper = 0;
    N = length(grid);
    
    %find the interval
    for i = 1:N-1
        if  x <= grid(i+1) 
            lower = i;
            upper = i+1;
            break;
        end
    end
    
    %get linear interpolant
    p = polyfit([grid(lower); grid(upper)], [u_h(lower); u_h(upper)], 1);
    uh = p(1)*x + p(2);
    uh_prime = p(1);
    