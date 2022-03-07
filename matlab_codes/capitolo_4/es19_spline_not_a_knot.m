x = linspace(-1,1,110);
fx = cos(pi*(x.^2)/2);

% Calcolo di una spline cubica not-a-knot con la function spline di Matlab
% su n + 1 ascisse equidistanti e su n + 1
% ascisse di Chebyshev nell'intervallo [-1, 1]
% e graficazione in forma semilogy del massimo 
% errore commesso.

max_errorSpline = zeros(20,1);
max_errorSplineCheby = zeros(20,1);

for n=4 : 100
    
    % Studio su  n + 1 ascisse equidistanti
    xi = linspace(-1,1,n+1);
    fi = cos(pi*(xi.^2)/2);
    s = spline(xi,fi,x);
    error = abs(fx - s);
    max_errorSpline(n) = max(error);
    
    % Studio su  n + 1 ascisse di Chebyshev
    xi_cheby = cheby(-1, 1, n+1);
    fi_cheby = cos(pi*(xi_cheby.^2)/2);
    s_cheby = spline(xi_cheby,fi_cheby,x);
    error2 = abs(fx - s_cheby);
    max_errorSplineCheby(n) = max(error2);
    
end
figure
semilogy(max_errorSpline);
figure
semilogy(max_errorSplineCheby);


function x = cheby(a,b,n)
    %
    % x=cheby(a,b,n) calcola le ascisse di Chebyshev per
    % il polinomio di interpolazione di grado n in [a,b]
    %
    if((a>=b) || n~=fix(n) || n<0)
        error("dati errati");
    end
    x = zeros(n,1);
    for i=1: n
        coseno = cos((2*(i-1)+1)*(pi)/(2*(n-1)+2));
        x(i) = (a+b)/2 + ((b-a)/2)*coseno;
    end
    x = flip(x);
    return
end
