x = linspace(-1,1, 110);
fx =  cos(pi*(x.^2)/2);

max_error_equidistanti = zeros(96,1);
max_error_cheby = zeros(96,1);

% Calcolo della spline cubica naturale 
% su n + 1 ascisse equidistanti e su n + 1
% ascisse di Chebyshev nell'intervallo [-1, 1]
% e graficazione in forma semilogy del massimo 
% errore commesso.

for n=4: 100
    
    % Studio su n + 1 ascisse equidistanti
    xi = linspace(-1,1,n+1);
    fi = cos(pi*(xi.^2)/2);
    spline = splinenat(xi,fi, x);
    error = abs(fx - spline);
    max_error_equidistanti(n) = max(error);
    
    % Studio su n + 1 ascisse di Chebyshev in [-1,1]
    xi_cheby = cheby(-1, 1, n+1);
    fi_cheby = cos(pi*(xi_cheby.^2)/2);
    spline2 = splinenat(xi_cheby,fi_cheby, x);
    error2 = abs(fx - spline2);
    max_error_cheby(n) = max(error2);
    
end

figure
semilogy(max_error_equidistanti);
figure
semilogy(max_error_cheby);



function spline = splinenat(xi, fi, xx)
    %
    %   spline = splinenat(xi, fi, xx) 
    %   Calcola la spline cubica naturale interpolante di una certa funzione
    %   
    %   La function prende in input:
    %   - xi: ascisse di interpolazione.
    %   - fi: valori assunti dalla funzione interpolata in corrispondenza delle
    %         ascisse di interpolazione.
    %   - xx: ascisse in cui avviene il calcolo della spline.
    %    
    %   La function restituisce in output:
    %    - spline: valori assunti dalla spline sulle ascisse xx.
    %

    mi= spline0(xi,fi);
    n= length(xi)-1;
    spline = zeros(size(xx));
    k= 1;
  
    for i=2:n+1 
        hi = xi(i) - xi(i-1);
        ri = fi(i-1) - (hi^2)/6 * mi(i-1);
        qi = (fi(i) - fi(i-1))/hi - hi/6 * (mi(i) - mi(i-1));
       
        while xx(k) < xi(i)
            spline(k) = (((xx(k) - xi(i-1))^3 * mi(i) + ...
                + (xi(i) - xx(k))^3 * mi(i-1))/(6 * hi) + ... 
                + qi * (xx(k) - xi(i-1)) + ri);
            k = k+1;
        end
    end
    return
end

function mi = spline0(xi, fi)
    %
    %   mi = spline0(xi, fi) Calcola il vettore degli mi per il calcolo di una
    %                        spline cubica naturale interpolante i punti (xi, fi).
    %
    m = length(xi);
    if m~=length(fi), error('dati errati'); end
    for i = 1:m-1
        if any( find(xi(i+1:m)==xi(i)) ), error('ascisse non distinte'), end
    end
    xi = xi(:); fi = fi(:);
    [xi,ind] = sort(xi); fi = fi(ind); % ordino le ascisse in modo crescente
    hi = diff(xi);
    n = m-1;
    df = diff(fi)./hi;
    hh = hi(1:n-1)+hi(2:n);
    rhs = 6*diff(df)./hh;
    phi = hi(1:n-1)./hh;
    csi = hi(2:n)./hh; % = 1-phi;
    d = 2*ones(n-1,1);
    phi = phi(2:n-1);
    csi = csi(1:n-2);
    mi = trisolve( phi, d, csi, rhs );
    mi = [0;mi;0];
    return
end

function x = trisolve( phi, d, csi, b )
    n = length(d);
    x = b;
    for i = 1:n-1
        phi(i) = phi(i)/d(i);
        d(i+1) = d(i+1) -phi(i)*csi(i);
        x(i+1) = x(i+1) -phi(i)*x(i);
    end
    
    x(n) = x(n)/d(n);
    for i = n-1:-1:1
        x(i) = (x(i)-csi(i)*x(i+1))/d(i);
    end
    return
end

function x = cheby(a,b,n)
    %
    %   x=cheby(a,b,n) calcola le ascisse di Chebyshev per
    %   il polinomio di interpolazione di grado n in [a,b]
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
