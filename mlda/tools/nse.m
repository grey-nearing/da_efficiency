function nse = nse(M,O)

% sgO = var(O);
muO = mean(O);
sqM = sum((M-O).^2);
sqO = sum((O-muO).^2);

nse = 1 - sqM/sqO;