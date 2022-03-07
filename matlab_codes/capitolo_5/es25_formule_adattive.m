f = @(x) 1./(1+(10.^2)*x.^2);

[integral, points, iterations] = adaptrap(f,-1,1,1e-6);
integral

[integral, points, iterations] = adapsim(f,-1,1,1e-6);
integral







function [I2,points, iterations] = adaptrap( f, a, b, tol, fa, fb )
%
%   function [I2,points, iterations] = adaptrap( f, a, b, tol )
%   Calcola un'approssimazione dell'integrale della funzione f tra gli
%   estremi di integrazione a,b con il metodo dei trapezi, 
%   concedendosi un errore di calcolo inferiore al valore di tol.
%
%   In input:
%   - f: funzione di cui calolare l'integrale
%   - a,b: estremi di integrazione
%   - tol: tolleranza utilizzata durante il calcolo
%   - fa, fb: funzione f calcolata negli estremi dei sottointervalli
%   calcolati
%
%   In output:
%   - I2 = valore dell'integrale
%   - points = punti in cui calcolo la funzione
%   - iterations = numero di chiamate ricorsive prima di ottenere
%   un'approssimazione soddisfacente
%

    global points
    global iterations
    delta = 0.5; % ampiezza minima intervalli
    if nargin<=4
        iterations = 0;
        fa = feval( f, a );
        fb = feval( f, b );
        if nargout==2
            points = [a fa; b fb];
        else
            points = [];
        end
    end
    iterations = iterations + 1;
    h = b-a;
    x1 = (a+b)/2;
    f1 = feval( f, x1 );
    if ~isempty(points)
        points = [points; [x1 f1]];
    end
    I1 = .5*h*( fa+fb );
    I2 = .5*( I1 + h*f1 );
    e = abs( I2-I1 )/3;
    if e>tol || abs(b-a)>delta
        I2 = adaptrap( f, a, x1, tol/2, fa, f1 ) +...
        adaptrap( f, x1, b, tol/2, f1, fb );
    end
    return
end

function [I2,points, iterations] = adapsim( f, a, b, tol, fa, f1, fb )
%
%   [I2,points, iterations] = adapsim( f, a, b, tol )
%
%   Calcola un'approssimazione dell'integrale della funzione f tra gli
%   estremi di integrazione a,b con il metodo di Simpson, 
%   concedendosi un errore di calcolo inferiore al valore di tol.
%
%   In input:
%   - f: funzione di cui calolare l'integrale
%   - a,b: estremi di integrazione
%   - tol: tolleranza utilizzata durante il calcolo
%   - fa, f1, fb: funzione f calcolata negli estremi dei sottointervalli
%   calcolati
%
%   In output:
%   - I2 = valore dell'integrale
%   - points = punti in cui calcolo la funzione
%   - iterations = numero di chiamate ricorsive prima di ottenere
%   un'approssimazione soddisfacente
%
    global points
    global iterations
    delta = 0.5; % ampiezza minima intervalli
    x1 = (a+b)/2;
    if nargin<=4
        iterations = 0;
        fa = feval( f, a );
        fb = feval( f, b );
        f1 = feval( f, x1 );
        if nargout==2
            points = [a fa;x1 f1; b fb];
        else
            points = [];
        end
    end
    iterations = iterations +1;
    h = (b-a)/6;
    x2 = (a+x1)/2;
    x3 = (x1+b)/2;
    f2 = feval( f, x2 );
    f3 = feval( f, x3 );
    if ~isempty(points)
        points = [points; [x2 f2; x3 f3]];
    end
    I1 = h*( fa+4*f1+fb );
    I2 = .5*h*( fa + 4*f2 + 2*f1 + 4*f3 +fb );
    e = abs( I2-I1 )/15;
    if e>tol || abs(b-a)>delta
        I2 = adapsim( f, a, x1, tol/2, fa, f2, f1 ) + ...
        adapsim( f, x1, b, tol/2, f1, f3, fb );
    end
    return
end