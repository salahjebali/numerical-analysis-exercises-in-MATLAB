function radice = secanti(func, iMax, x0, x, tol)
%
%   radice = secanti(func, iMax, x0, x, tol)
%   calcola la radice di una funzione con il metodo delle Secanti.
%   
%   La function prende in input: 
%   -func: funzione di cui calcolare la radice.
%   -iMax: numero di iterazioni massimo.
%   -x0,x: approssimazioni iniziali.
%   -tol: tolleranza per criterio di arresto richiesto
%
%   La function resituisce in output:
%   - [radice, iterations] = [radice, numero di iterazioni per ottenerla]
%

    fx0=feval(func,x0);  
    fx = feval(func,x);
    x1 = (fx*x0 - fx0*x)/(fx - fx0);
    iterations = 1; 
    while (iterations < iMax) && (abs(x1-x) > ( tol * (1 + abs(x)) ))
        x0 = x;
        x = x1;
        iterations= iterations+1;
        fx0 = fx;
        fx = feval(func,x);
        x1 = (fx*x0 - fx0*x)/(fx - fx0);
    end
    if abs(x1-x) > (tol * (1 + abs(x)))
         disp('Il metodo non converge');
    end
    radice = x1;
    return
end
