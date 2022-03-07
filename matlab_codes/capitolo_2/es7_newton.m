f=@func;
f1=@der;
tol= 1e-3; % 1e-6, 1e-9, 1e-12
x0 = 1;
[x,i] = newton(f, f1, 100, x0, tol);

function y = func(x)
    y=(x^2)*tan(x);
end

function y = der(x)
    y = (x^2)/(cos(x))^2 + 2*x*tan(x);
end

function [x,iterations] = newton(func, f1, iMax, x0, tol)
%
%   [x,iterations] = newton(func, f1, iMax, x0, tol)
%   calcola la radice di una funzione con il metodo di Newton.
%   
%   La function prende in input: 
%
%   -func: funzione di cui calcolare la radice.
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
    fx0=feval(func,x0);
    f1x=feval(f1,x0);
    x=x0-(fx0/f1x);    
    while (iterations < iMax) && (abs(x-x0) > ( tol * (1 + abs(x0)) ) )
        iterations= iterations+1;
        x0=x;
        fx0= feval(func,x0);
        f1x= feval(f1,x0);
        x=x0-(fx0/f1x);
    end
    if abs(x-x0) > (tol * (1 + abs(x0)))
         disp('il metodo non converge');
    end
    return
end