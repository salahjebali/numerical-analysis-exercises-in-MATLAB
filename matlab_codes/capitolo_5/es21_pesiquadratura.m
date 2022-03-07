
c = calcolaPesi(7)



function c = calcolaPesi(n)
%
%   c = calcolaPesi(n) calcola i pesi della quadratura
%   della formula di Newton-Cotes di grado n.
%
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




