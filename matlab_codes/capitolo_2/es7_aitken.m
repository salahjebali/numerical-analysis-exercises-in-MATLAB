f=@func;
f1=@der;
tol= 1e-3; % 1e-6, 1e-9, 1e-12
x0 = 1;
[x,i] = Aitken(f, f1, 100, x0, tol);

function y = func(x)
    y=(x^2)*tan(x);
end

function y = der(x)
    y = (x^2)/(cos(x))^2 + 2*x*tan(x);
end

function [x,iterations] = Aitken(f, f1, iMax, x0, tol)
%
%   [x,iterations] = Aitken(f, f1, iMax, x0, tol)
%   calcola la radice di una funzione con il metodo di Aitken.
%   
%   La function prende in input: 
%
%   -f: funzione di cui calcolare la radice.
%   -f1: derivata prima della funzione.
%   -iMax: numero di iterazioni massimo.
%   -x0: approssimazione iniziale.
%   -tol: tolleranza per criterio di arresto richiesto
%
%   La function restituisce in output: 
%
%   - [x, iterations] = [radice, numero di iterazioni per ottenerla]
%
    iterations = 0;
    fx = feval(f,x0);
    if fx==0
        x=x0;
        return 
    end
    f1x = feval(f1,x0); 
    x= x0-fx/f1x;
    vai = 1;
    while (iterations<iMax) && vai
        iterations = iterations + 1;
        x0 = x;
        fx = feval(f, x0);
        f1x = feval(f1, x0);
        x1 = x0 - fx/f1x;
        fx = feval(f,x1);
        f1x = feval(f1,x1);
        x = x1 - fx/f1x;
        x = (x*x0-x1^2)/(x-2*x1+x0);
        vai = (abs(x-x0) > ( tol * (1 + abs(x0))));
    end
    if vai, disp('il metodo non converge'), end
end
