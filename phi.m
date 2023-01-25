function y =  phi(h, i , x)
    grid = cumsum(h);  
    if i == 1
        if 0 < x && x <= grid(i)
            y = x / h(i);
        elseif grid(1) < x && x <= grid(i+1)
            y = (grid(i+1)-x)/h(i+1);
        else
            y = 0;
        end
    end
    
    if i > 1
        if grid(i-1) < x && x <= grid(i)
            y = (x-grid(i-1))/h(i);
        elseif grid(i) < x && x <= grid(i+1)
            y = (grid(i+1)-x)/h(i+1);
        else
            y = 0;
        end
    end
    

    