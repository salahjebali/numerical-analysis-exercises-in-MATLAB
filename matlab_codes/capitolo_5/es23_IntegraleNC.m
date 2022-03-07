I_f = log(cos(1)/cos(1.1));

% Calcolo dell'integrale di una funzione mediante formule di Newton-Cotes 
% di grado n=1..9 
for n = 1:9
    I_nc = NewtonCotes(@tan, n, -1, 1.1);
    fprintf('n = %d\t', n);
    fprintf('In(f) = %.3d\t', I_nc);
    fprintf('errore = %.3d\n', abs(I_f - I_nc));  
end

function I_nc = NewtonCotes(func, n, a, b)
%
%   I_nc = NewtonCotes(func, n, a, b) calcola l'approssimazione 
%   dell' integrale definito tra "a" e "b" attraverso la formula
%   di quadratura di Newton-Cotes.
%   
%   La function prende in input: 
%   - func = funzione integranda.
%   - n = grado della fomrula di Newton-Cotes.
%   - a,b = estremi di integrazione.
%   
%   La function restituisce in output:
%   - I_nc = Approssimazione integrale ottenuta da Newton-Cotes
%
    h = (b - a) / n;
    x = a + (0:n) * h;
    fx = func(x);
    c = calcolaPesi(n);
    I_nc = 0;
    for i = 1:n+1
        I_nc = I_nc + c(i) * fx(i);
    end
    I_nc = I_nc * h;
end

function c = calcolaPesi(n)
%
%   c = calcolaPesi(n) calcola i pesi della quadratura
%   della formula di Newton-Cotes di grado n.
%   La function prende in input: 
%   - n = grado della formula di quadratura Newton-Cotes.
%
%   La function restituisce in output: 
%   - c = vettore dei pesi.
%
    c = zeros(1,n+1);
    for i=1:n+1
        c(i) = integral(@calcolaLin, 0, n);
    end
    
    function Lin = calcolaLin(t) 
        indexi = i-1;
        ind = [0:indexi-1 indexi+1:n];
        Lin = 1/prod( indexi-ind );
        for j = 1:n
            Lin = Lin.*( t - ind(j));
        end
    end
end