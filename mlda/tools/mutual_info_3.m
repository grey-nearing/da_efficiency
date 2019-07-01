function [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,Bx,By,Bz)

N = length(X);
assert(N==length(Y))
assert(N==length(Z))

[Pxyz,Pxy,Pxz,Pyz,Px,Py,Pz] = hist3(X,Y,Z,Bx,By,Bz);

Px = Px(:);
Py = Py(:);
Pz = Pz(:);
Pyz = Pyz(:);
Pxz = Pxz(:);
Pxy = Pxy(:);
Pxyz = Pxyz(:);

if abs(sum(Px)-1)>1/N^2;   error('Px does not sum to 1');   end;
if abs(sum(Py)-1)>1/N^2;   error('Py does not sum to 1');   end;
if abs(sum(Pz)-1)>1/N^2;   error('Pz does not sum to 1');   end;
if abs(sum(Pxy)-1)>1/N^2;  error('Pxy does not sum to 1');  end;
if abs(sum(Pxz)-1)>1/N^2;  error('Pxz does not sum to 1');  end;
if abs(sum(Pyz)-1)>1/N^2;  error('Pyz does not sum to 1');  end;
if abs(sum(Pxyz)-1)>1/N^2; error('Pxyz does not sum to 1'); end;

Hx   = -Px(Px>0)'     * log(Px(Px>0));
Hy   = -Py(Py>0)'     * log(Py(Py>0));
Hz   = -Pz(Pz>0)'     * log(Pz(Pz>0));
Hxy  = -Pxy(Pxy>0)'   * log(Pxy(Pxy>0));
Hxz  = -Pxz(Pxz>0)'   * log(Pxz(Pxz>0));
Hyz  = -Pyz(Pyz>0)'   * log(Pyz(Pyz>0));
Hxyz = -Pxyz(Pxyz>0)' * log(Pxyz(Pxyz>0));

Ixy  = Hx+Hy-Hxy;
Ixz  = Hx+Hz-Hxz;
Iyz  = Hy+Hz-Hyz;
Ixyz = -(Hx+Hy+Hz-Hxyz-Ixy-Ixz-Iyz);









