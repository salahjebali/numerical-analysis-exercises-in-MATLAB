function x = qrsolve(QR, b)
%
%   x = qrsolve(QR, b) 
%   Data in ingresso la matrice fattorizzata QR e il
%   termine noto b, risolve il sistema Ax = b nel senso dei minimi quadrati.
%

    [m,n] = size(QR);
    k = length(b);
    
    if k ~= m || n > m
        error('Dati Inconsistenti');
    end
    
    x = b(:);
    
    for i= 1 : n
        v = [1; QR(i+1:m, i)];
        beta = 2/(v' * v);
        x(i:m) = x(i:m) - (beta * (v' * x(i:m)))*v;
    end
    
    x = x(1:n);
    
    for i=n:-1:1
        if (QR(i,i) == 0)
            error('Sistema irrisolvibile');
        end
        
        x(i) = x(i)/QR(i,i);
        if i>1
            x(1:i-1) = x(1:i-1) - x(i) * QR(1:i-1, i);
        end
    end
    return
end