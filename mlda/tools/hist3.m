function [Pxyz,Pxy,Pxz,Pyz,Px,Py,Pz] = hist3(X,Y,Z,Bx,By,Bz)

N = length(X);
assert(N==length(Y))
assert(N==length(Z))

Mx = length(Bx);
My = length(By);
Mz = length(Bz);

Pxyz = zeros(Mx,My,Mz);
for n = 1:N
  x = X(n);
  y = Y(n);
  z = Z(n);
  y = find(y>=By,1,'last');
  if isempty(y); keyboard; end
  x = find(x>=Bx,1,'last');
  if isempty(x); keyboard; end
  z = find(z>=Bz,1,'last');
  if isempty(z); keyboard; end
  Pxyz(x,y,z) = Pxyz(x,y,z) + 1;
end
%Pxyz(Pxyz==0)=1/1000;
Pxyz = Pxyz/sum(Pxyz(:));

Pxy = squeeze(sum(Pxyz,3));
Pxz = squeeze(sum(Pxyz,2));
Pyz = squeeze(sum(Pxyz,1));
Px  = squeeze(sum(Pxy,2));
Py  = squeeze(sum(Pxy,1));
Pz  = squeeze(sum(Pxz,1));








