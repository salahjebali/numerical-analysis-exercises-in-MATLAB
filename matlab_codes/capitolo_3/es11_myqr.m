function QR = myqr(A)
%
%   QR = myqr(A)
%   Data in ingresso una matrice A di numeri reali, la cui dimensione risulta
%   n x m con m >= n = rank(A), restituisca una matrice QR che contenga 
%   l'informazione sui fattori Q ed R della fattorizzazione QR di A.
%

    [m,n] = size(A);
    
    % Controlli di integrità 
    if n > m
        error('Dimensioni matrice non corrette');
    end
    
    QR = A;
    
    for i=1 : n
        alpha = norm(QR(i:m, i));
        
        if (alpha == 0)
            error('La matrice non è di rango massimo')
        end
        
        if (QR(i,i) >= 0)
            alpha = -alpha;
        end
        
        v1 = QR(i,i) - alpha;
        QR(i,i) = alpha;
        QR(i+1:m, i) = QR(i+1:m, i) / v1;
        beta = -v1/alpha;
        v = [1; QR(i+1:m, i)];
        QR(i:m,i+1:n) = QR(i:m, i+1:n) - (beta*v)*(v'*QR(i:m,i+1:n));
    end
    return
end