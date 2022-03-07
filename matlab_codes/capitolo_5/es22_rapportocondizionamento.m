rapporto_k = rapporto_condizionamento(-5,5);
semilogy(rapporto_k);

function rapporto_k = rapporto_condizionamento(a,b)
    %
    %   rapporto_k = rapporto_condizionamento(a,b) calcola
    %   il rapporto kn/k, ovvero il rapporto tra il numero
    %   di condizionamento della formula di quadratura di 
    %   Newton-Cotes e quello del calcolo dell'integrale 
    %   definito tra "a" e "b".
    %
    %   La function prende in input:
    %   - a,b = estremi di integrazione.
    %
    %   La function restituisce in output:
    %   - rapporto_k = vettore contenente valori di kn/k.
    %
    
    k = b-a;

    for n = 1:50
        c  = calcolaPesi(n);
        kn = (b-a)*(sum(abs(c))/n);
        rapporto_k(n) = kn/k;
    end

    return
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
    return
end

