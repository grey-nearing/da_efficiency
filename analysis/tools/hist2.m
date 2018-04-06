function [Pxy,Px,Py] = hist2(X,Y,Bx,By)

N = length(X);
assert(N==length(Y))

Mx = length(Bx);
My = length(By);

Pxy = zeros(Mx,My);
for n = 1:N
  x = X(n);
  y = Y(n);
  y = find(y>=By,1,'last');
  if isempty(y); keyboard; end
  x = find(x>=Bx,1,'last');
  if isempty(x); keyboard; end
  Pxy(x,y) = Pxy(x,y) + 1;
end
%Pxy(Pxy==0)=1/1000;
Pxy = Pxy/sum(Pxy(:));

Px  = squeeze(sum(Pxy,2));
Py  = squeeze(sum(Pxy,1));








