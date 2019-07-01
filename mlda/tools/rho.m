function r2 = rho(X,Y)

c = corrcoef(X,Y);
r2 = c(1,2);
