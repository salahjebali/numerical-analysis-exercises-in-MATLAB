function radice = bisezione(func, a, b, tol)
%   radice = bisezione(func, a, b, tol) 
%   calcola la radice di una funzione con il metodo di bisezione. 
%   La function prende in input:
%
%   -func: funzione di cui calcolare la radice.
%   -a, b: estremi intervallo di confidenza iniziale.
%   -tol: tolleranza per criterio di arresto richiesto.
%
%   La function resituisce in output:
%   - [radice, iterations] = [radice, numero di iterazioni per ottenerla]
%
    fa=feval(func,a);
    imax= ceil( log2(b-a) - log2(tol) );
    x2= (a+b)/2;
    for iterations=1:imax
        x = x2;
        fx = feval(func,x);
        if fa*fx > 0
            a=x;
            fa=fx;
        else
            b=x;
        end
        x2= (a+b)/2;
        if abs(x2 - x) <= (tol * (1 + abs(x)))
            break;
        end
    end
    radice = x2;
    return
end
