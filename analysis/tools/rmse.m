function r = rmse(X,Y,norm);

[N,D] = size(X);
assert(N == size(Y,1));
assert(D == size(Y,2));

X = X-mean(X);
Y = Y-mean(Y);

for d = 1:D
 r(d) = sqrt(mean((X(:,d)-Y(:,d)).^2));
 if     nargin == 3 && norm == 1
  r(d) = r(d) / sqrt(mean(X(:,d).^2));
 elseif nargin == 3 && norm == 2
  r(d) = r(d) / sqrt(mean(Y(:,d).^2));
 end
end


