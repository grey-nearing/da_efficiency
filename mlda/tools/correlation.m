function r2 = correlation(X,Y)

c = corrcoef(X,Y);
r2 = c(1,2);
