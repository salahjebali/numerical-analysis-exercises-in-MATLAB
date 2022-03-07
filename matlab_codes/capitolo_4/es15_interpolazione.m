x = linspace(-1,1,100); 
fx = cos(pi*(x.^2)/2);

% Calcolo del polinomio in forma di Lagrange
% su n + 1 ascisse equidistanti e su n + 1
% ascisse di Chebyshev nell'intervallo [-1, 1]
% e graficazione in forma semilogy del massimo 
% errore commesso.


max_error_equidistanti = zeros(40,1);
max_error_cheby = zeros(40,1);

for n = 1:40  
    
    % Studio su  n + 1 ascisse equidistanti
    xn = linspace(-1,1,n+1); 
    fn = cos(pi*(xn.^2)/2);
    yl = lagrange(xn, fn, x);
    error = abs(fx-yl);
    max_error_equidistanti(n) = max(error);
    
    % Studio su n + 1 ascisse di Chebyshev in [-1,1]
    xn_cheby = cheby(-1,1,n+1);
    fn_cheby = cos(pi*(xn_cheby.^2)/2);
    y_cheby = lagrange(xn_cheby, fn_cheby, x);
    error = abs(fx-y_cheby);
    max_error_cheby(n) = max(error);
    
end

figure
semilogy(max_error_equidistanti);
figure
semilogy(max_error_cheby);

function y = lagrange( xi, fi, x )
    %
    % y = lagrange( xi, fi, x ) calcola il polinomio interpolante sulle
    % coppie (xi,fi), in forma di Lagrange, valutato nei punti del vettore x.
    %

    n = length(xi); % il grado del polinomio interpolante `e n-1
    if length(fi)~=n, error('dati incosistenti');
    else
    for i = 1:n-1
        if any(find(xi(i+1:n)==xi(i))),error('ascisse non distinte'),end
        end
    end
    y = zeros(size(x));
    for i = 1:n
        y = y + fi(i)*Li(i,xi,x);
    end
    return
end

function Lx = Li(i,xi,x)
    %
    % Lx = Li(i,xi,x) calcola l'i-esimo elemento 
    % della base di Lagrange.
    %
    c = xi;
    ci = c(i);
    c(i) = []; 
    Lx = ones(size(x));

    for k = 1:length(c)
        Lx = Lx.*(x-c(k));
    end

    Lx = Lx/prod(ci-c);

    return
end

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
