n = 10^4;
xi = -1+2*(0:n)/ n;
fi =  cos(pi.*(xi)/2);
% Calcolo fi perturbata a partire dal valore reale di fi
fipert = fi + 10^-3 * rand(size(xi));

errore = zeros(20);
for m = 1:20
    % Calcolo in ordine descrescente di potenza i coefficienti del
    % polinomio interpolante
    a = polyfit(xi, fipert, m);
    % Calcolo il polinomio sulle ascisse di interpolazione
    poli = polyval(a, xi);
    % Calcolo la norma euclidea della differenza tra la funzione calcolata
    % nelle ascisse di interpolazione ed il polinomio trovato
    errore(m) = norm(fi - poli);
end
semilogy(errore);

