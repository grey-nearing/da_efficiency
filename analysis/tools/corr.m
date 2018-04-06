function c = corr(X,Y);

[N,D] = size(X);
assert(N == size(Y,1));
assert(D == size(Y,2));

for d = 1:D
 cc = corrcoef(X(:,d),Y(:,d));
 c(d) = cc(2);
end


