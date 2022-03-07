f = @(x) tan(x);
a = -1;
b = 1.1;

for i=1:9
    T = TrapeComp(f, a, b, i);
end

for i=2:2:8
    S = SimpComp(f, a, b, i)
end

function I=TrapeComp(fun,a,b,n)
    %
    %   I = TrapeComp(fun,a,b,n) calcola l'approssimazione dell'integrale 
    %   definito di f(x) con estremi "a" e "b", mediante la formula
    %   composita dei trapezi su n+1 ascisse equidistanti.
    %
    %   La function prende in input:
    %   - fun = stringa con nome della function che implementa la funzione
    %   integranda che accetta input vettoriali.
    %   - a,b = estremi di integrazione.
    %   - n = numero di sottointervalli.
    %   
    %   La function restituisce in output:
    %   - I = valore dell'approssimazione dell'integrale.
    %
    if a==b
        I = 0;
    elseif(n < 1 || n ~= fix(n))
        error("numero ascisse non valido");
    else
        h = (b-a)/n;
        x = linspace(a,b,n+1);
        f = feval(fun,x);
        I = h*(f(1)/2 + sum(f(2:n)) + f(n+1)/2);
    end
    return
end


function I = SimpComp(fun,a,b,n)
    %
    %   I = SimpComp(fun,a,b,n) calcola l'approssimazione dell'integrale 
    %   definito di f(x) con estremi "a" e "b", mediante la formula
    %   composita di Simpson su n+1 ascisse equidistanti.
    %
    %   La function prende in input:
    %   - fun = stringa con nome della function che implementa la funzione
    %   integranda che accetta input vettoriali.
    %   - a,b = estremi di integrazione.
    %   - n = numero di sottointervalli.
    %   
    %   La function restituisce in output:
    %   - I = valore dell'approssimazione dell'integrale.
    %
    if a==b
        I = 0;
    elseif(n < 2 || n/2 ~= fix(n/2))
        error("numero ascisse non valido");
    else
        h = (b-a)/n;
        x = linspace(a,b,n+1);
        f = feval(fun,x);
        I = (h/3)*(f(1) + f(n+1) +...
            4*sum(f(2:2:n)) + ...
            2*sum(f(3:2:n-1)));
    end
    return
end