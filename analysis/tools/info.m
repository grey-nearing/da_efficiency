function i = info(X,Y,Bx,By,norm)

N = length(X);
assert(N == size(Y,1));

[N,D] = size(X);
assert(N == size(Y,1));
assert(D == size(Y,2));

for d = 1:D
 [i(d),hx,hy] = mutual_info(X(:,d),Y(:,d),Bx,By);
 if     nargin == 5 && norm == 1
  i(d) = i(d) / hx;
 elseif nargin == 5 && norm == 2
  i(d) = i(d) / hy;
 end
end


