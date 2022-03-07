function radice = corde(func, der, iMax, x0, tol)
%
%   radice = corde(func, der, iMax, x0, tol)
%   calcola la radice di una funzione con il metodo delle Corde.
%   
%   La function prende in input: 
%   -func: funzione di cui calcolare la radice.
%   -der: derivata prima della funzione.
%   -iMax: numero di iterazioni massimo.
%   -x0: approssimazione iniziale.
%   -tol: tolleranza per criterio di arresto richiesto.
%
%   La function resituisce in output:
%   - [radice, iterations] = [radice, numero di iterazioni per ottenerla]
%
    iterations = 0;
    fx0=feval(func,x0);
    derx0=feval(der,x0);
    x=x0-(fx0/derx0);    
    while (iterations < iMax) && (abs(x-x0) > ( tol * (1 + abs(x0)) ) )
        iterations= iterations+1;
        x0=x;
        fx0= feval(func,x0);
        x=x0-(fx0/derx0);
    end
    if abs(x-x0) > (tol * (1 + abs(x0)))
         disp('Il metodo non converge');
    end
    radice = x;
    return
end
