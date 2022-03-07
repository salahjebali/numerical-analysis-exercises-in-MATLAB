f=@func;
f1=@der;
tol= 1e-3; % 1e-6, 1e-9, 1e-12
x0 = 1;
m = 3;
[x,i] = NewtonModificato(f, f1, m, 100, x0, tol);

function y = func(x)
    y=(x^2)*tan(x);
end

function y = der(x)
    y = (x^2)/(cos(x))^2 + 2*x*tan(x);
end

function [x,iterations] = NewtonModificato(func, f1, m, iMax, x0, tol)
%
%   [x, iterations] = NewtonModificato(func, f1, m, iMax, x0, tol)
%   calcola la radice di una funzione con il metodo Newton Modificato.
%   La function prende in input:
%   
%   - func = function di cui calcolare la radice.
%   - f1 = derivata prima della funzione.
%   - m = molteplicità della radice.
%   - iMax = numero di iterazioni massime.
%   - x0 = approssimazione radice iniziale.
%   - tol = tolleranza.
%  
%   La function resituisce in output:
%   - [x, iterations] = [radice, numero di iterazioni per ottenerla]
%
    iterations = 0;
    fx0=feval(func,x0);
    f1x=feval(f1,x0);
    x=x0-m*(fx0/f1x);    
    while (iterations < iMax) && (abs(x-x0) > ( tol * (1 + abs(x0)) ) )
        iterations= iterations+1;
        x0=x;
        fx0= feval(func,x0);
        f1x= feval(f1,x0);
        x=x0-m*(fx0/f1x);
    end
    if abs(x-x0) > (tol * (1 + abs(x0)))
         disp('il metodo non converge');
    end
end