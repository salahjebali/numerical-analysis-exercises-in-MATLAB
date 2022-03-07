x = linspace(-1,1,100);
fx = cos(pi*(x.^2)/2);

% Calcolo del polinomio di Hermite
% su n + 1 ascisse equidistanti e su n + 1
% ascisse di Chebyshev nell'intervallo [-1, 1]
% e graficazione in forma semilogy del massimo 
% errore commesso.

max_errorHermite = zeros(20,1);
max_errorHermiteCheby = zeros(20,1);
for n=1 : 20
    
    % Studio su  n + 1 ascisse equidistanti
    xi = linspace(-1,1,n+1);
    fi = cos(pi*(xi.^2)/2);
    f1i = -(pi*xi).*sin(pi*(xi.^2)/2); 
    yh = hermite(xi,fi,f1i,x);
    error = abs(fx - yh);
    max_errorHermite(n) = max(error);
    
    % Studio su n + 1 ascisse di Chebyshev in [-1,1]
    xi_cheby = cheby(-1, 1, n+1);
    fi_cheby = cos(pi*(xi_cheby.^2)/2);
    f1i_cheby = -(pi*xi_cheby).*sin(pi*(xi_cheby.^2)/2);
    yh2 = hermite(xi_cheby,fi_cheby,f1i_cheby,x);
    error2 = abs(fx - yh2);
    max_errorHermiteCheby(n) = max(error2);
    
end
figure
semilogy(max_errorHermite);
figure
semilogy(max_errorHermiteCheby);


function [y,df] = hermite( xi, fi, f1i, xx )
% [y,df] = hermite( xi, fi, f1i, xx ) Calcola il valore del polinomio
% interpolante di Hermite sulle
% ascisse xi. I vettori fi e f1i
% contengono i corrispondenti valori della f e della sua ferivata prima.
% Se specificato, df contiene il vettore delle differenze divise.

%
% controlli sui dati di ingresso
%
    m = length(xi);
    if m~=length(fi) || m~=length(f1i), error('dati inconsistenti'), end
    for i = 1:m-1
        if any( find(xi(i+1:m)==xi(i)) ), error('ascisse non distinte'), end
    end

    n = 2*m-1; % grado del polinomio interpolante
    x = zeros(n+1,1);
    df = x;
    x(1:2:n) = xi(:);
    x(2:2:n+1) = xi(:);
    df(1:2:n) = fi(:);
    df(2:2:n+1) = f1i(:);

    for i = n:-2:3 % seconda colonna della tabella
        df(i) = ( df(i)-df(i-2) )/( x(i)-x(i-1) );
    end
    for i = 2:n % colonne successive della tabella
        for j = n+1:-1:i+1
            df(j) = ( df(j)-df(j-1) )/( x(j)-x(j-i) );
        end
    end

    % calcolo il polinomio interpolante nelle ascisse prescritte
    %
    y = df(n+1)*ones(size(xx));
    for k = 0:n-1
        y = y.*( xx-x(n-k) ) +df(n-k);
    end
    return
end


function x = cheby(a,b,n)
%   x=cheby(a,b,n) calcola le ascisse di Chebyshev per
%   il polinomio di interpolazione di grado n in [a,b]
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
