function radice = radn(x, n)
%
%   radice = radn(x,n) calcola la radice n-esima di un numero positivo x.
%   Il problema del calcolo della radice n-esima di x si riconduce a cerca
%   la radice della funzione f(x) = x^(1/n) - radice. In questo caso è
%   stato adottato il metodo di Newton poiché il più efficiente tra quelli
%   studiati.
%   
%   La function prende in input: 
%   - x: numero positivo di cui si vuole calcolare la radice.
%   - n: numero intero che indica il grado della radice.
%
%   In output viene restituita la radice calcolata.
%

x_start = x;    
x0 = x;
fx = f(x0,n,x_start);
derx = der(x0,n);
x = x0 - fx/derx;
while (abs(x-x0)>eps(x))
    x0 = x;
    fx = f(x0,n,x_start);
    derx = der(x0,n);
    x = x0- fx/derx;
end
radice = x;    
end

function y = f(x, n, x_start)
y = x^n - x_start;
end

function y = der(x, n)
y = n * x ^ (n-1);
end