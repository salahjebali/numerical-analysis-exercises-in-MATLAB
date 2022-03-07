function [LU,p] = palu(A)
%
%   [LU,p] = palu(A) 
%   Calcola la fattorizzazione LU di una matrice quadrata A con a tecnica
%   del pivoting parziale, restituendo il vettore delle permutazioni
%   eseguite p.
%
%   La function prende in input:
%   - A: matrice quadrata non singolare, se la matrice passata non ha
%   questi requisiti viene rifiutata dalla funzione con messaggio di errore
%
%  In output vengono restituiti:
%  - LU: matrice contenente l'informazione sui fattori L ed U in cui è
%  stata scomposta la matrice di partenza
%  - p: vettore contenente le permutazioni effettuate sulla matrice durante
%  la fattorizzazione con tecnica del pivoting parziale

    LU = A;  
    [n,m] = size(A);
    if  (n~=m)
        error("Matrice non quadrata.");
    end
    for i=1:n-1
        if (LU(i,i) == 0)
            error("Matrice singolare, quindi non fattorizzabile LU.");
        end
    end


    p = 1:n;
    for i=1: n-1
        [mi,ki] = max(abs(LU(i:n,i)));
        if (mi==0)
            error("Matrice singolare, quindi non fattorizzabile LU.");
        end
        ki = ki + i - 1;
        if (ki>i)  
            p([i,ki]) = p([ki,i]);
            LU([i,ki],:) = LU([ki,i],:); 
        end
        LU(i+1:n,i) = LU(i+1:n,i)/LU(i,i);
        LU(i+1:n,i+1:n) = LU(i+1:n,i+1:n) - LU(i+1:n,i)*LU(i,i+1:n);
    end 
    return
end