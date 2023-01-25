hs = [1/10;1/20;1/40;1/80;1/160];
l2_norm = zeros(length(hs),1);
h1_seminorm = zeros(length(hs),1);

for i = 1: length(hs)
    h = hs(i)*(ones(1/hs(i) , 1));
    grid = cumsum(h);
    u_h = fin_elt_1d(h);
    N = length(h);
    
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
        
        %accumulate values for gaussian quadrature
        x1 = (-1/sqrt(3))*(b-a)/2 + (a+b)/2;
        x2 = (1/sqrt(3))*(b-a)/2 + (a+b)/2;
        exact_1 = exact( x1 );
        [approx_1, approx_1_slope] = u_heval(u_h, grid(1:end-1), x1 );
        exact_2 = exact( x2 );
        [approx_2, approx_2_slope] = u_heval(u_h, grid(1:end-1), x2 );
        l2_norm(i) = l2_norm(i)+ ((b-a)/2)*abs(exact_1-approx_1)^2 + ((b-a)/2)*abs(exact_2-approx_2)^2;
        
        exact_deriv_1 = derivative( x1 );
        exact_deriv_2 = derivative( x2 );
        h1_seminorm(i) = h1_seminorm(i) + ((b-a)/2)*abs(exact_deriv_1-approx_1_slope)^2 + ((b-a)/2)*abs(exact_deriv_2-approx_2_slope)^2;
        
    end
    l2_norm(i) = sqrt(l2_norm(i));
    h1_seminorm(i) = sqrt(h1_seminorm(i));
end

t = tiledlayout(1,2);
ax1 = nexttile; 
plot(log(hs), log(l2_norm), '-*', 'LineWidth', 2);
hold on
X = linspace(log(hs(1)),log(hs(end)), 100);
plot(X, 2*X-1, 'Linewidth', 2); 
legend('||u - u_{h}||_{L^{2}}', 'y = 2x-1', 'Location', 'best'); 
title('L^{2} Norm');
xlabel('h (log scale)');
ylabel('err(h) (log scale)')
hold off

ax2 = nexttile;
plot(log(hs), log(h1_seminorm), '-*', 'LineWidth', 2);
hold on
X = linspace(log(hs(1)),log(hs(end)), 100);
plot(X, X, 'Linewidth', 2); 
legend('|u - u_{h}|_{H^{1}}', 'y = x', 'Location', 'southeast'); 
title('H^{1} Seminorm (Energy Norm)');
xlabel('h (log scale)');
hold off
linkaxes([ax1, ax2]);