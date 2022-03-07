function x = lusolve(LU, p, b)
%
%   x = lusolve(LU, p, b)
%   Calcola la soluzione x del sistema lineare Ax=b in cui la matrice A è
%   stata fattorizzata mediante l'utilizzo della fattorizzazione LU con
%   tecnica di pivoting parziale.
%
%   La function prende in input:
%   - LU: matrice contenente l'informazione sui fattori L ed U in cui è
%  stata scomposta la matrice dei coefficienti del sistema lineare Ax=b
%   - p: vettore contenente le permutazioni effettuate sulle righe della 
%   matrice fattorizzata 
%   - b: termine noto del sistema lineare Ax=b
%
%  In output viene restituita:
%  - x: vettore contenente la soluzione del sistema lineare Ax=b
%

    [m,n] = size(LU);
    if ((m~=n) || (n~=length(b)) || (n~=length(p)))
        error('Dati inconsistenti.');
    elseif (min(abs(diag(LU)))== 0)
            error('Elemento sulla diagonale = 0.');
    end

    % Permuto gli elementi di b
    x = b(:);
    x = x(p);
    % Risolvo L
    for i=1: n-1 
        x(i+1:n) = x(i+1:n)-LU(i+1:n, i) * x(i);
    end
    % Risolvo U
    for i=n: -1: 1
        x(i) = x(i)/LU(i,i);
        x(1:i-1) = x(1:i-1) - LU(1:i-1, i) * x(i);
    end
    return
end