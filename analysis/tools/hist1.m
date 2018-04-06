function Px = hist1(X,Bx)

N = length(X);
Mx = length(Bx);

Px = zeros(Mx,1);
for n = 1:N
  x = X(n);
  x = find(x>=Bx,1,'last');
  if isempty(x); keyboard; end
  Px(x) = Px(x) + 1;
end
Px = Px/N;









