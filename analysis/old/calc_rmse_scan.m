clear all
close all
clc
addpath('tools');

load('data/all_data_scan.mat');

% ------ ALL SITES -----------------------
Xo = outdata.Mopen(:);
Xe = outdata.Moens(:);
Xa = outdata.Menks(:);
Xr = outdata.Mreg(:);
Yy = outdata.CDF(:);
Zz = outdata.Val(:);

assert(isempty(find(Xo<0))); 
assert(isempty(find(Xe<0))); 
assert(isempty(find(Xa<0))); 
assert(isempty(find(Xr<0))); 
assert(isempty(find(Zz<0))); 

I = find(isnan(Xo)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];  
I = find(isnan(Xe)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
I = find(isnan(Xa)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
I = find(isnan(Xr)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
I = find(isnan(Zz)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
I = find(isnan(Yy)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];

N = length(Zz)
rmse_all = [rmse(Xo,Zz),rmse(Xe,Zz),rmse(Xa,Zz),rmse(Xr,Zz),rmse(Yy,Zz)];
corr_all = [corr(Xo,Zz),corr(Xe,Zz),corr(Xa,Zz),corr(Xr,Zz),corr(Yy,Zz)];

% ------ BY SITE --------------------------
for s= 1:length(outdata.Ndata)

 Xo = outdata.Mopen(:,s);
 Xe = outdata.Moens(:,s);
 Xa = outdata.Menks(:,s);
 Xr = outdata.Mreg(:,s);
 Yy = outdata.CDF(:,s);
 Zz = outdata.Val(:,s);

 I = find(isnan(Xo)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];  
 I = find(isnan(Xe)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
 I = find(isnan(Xa)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
 I = find(isnan(Xr)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
 I = find(isnan(Zz)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];
 I = find(isnan(Yy)); Xo(I) = []; Xe(I) = []; Xa(I) = []; Zz(I) = []; Yy(I) = []; Xr(I) = [];

 N = length(Zz)
 rmseS(s,:) = [rmse(Xo,Zz),rmse(Xe,Zz),rmse(Xa,Zz),rmse(Xr,Zz),rmse(Yy,Zz)]
 corrS(s,:) = [corr(Xo,Zz),corr(Xe,Zz),corr(Xa,Zz),corr(Xr,Zz),corr(Yy,Zz)];

end

rmse_all
nanmean(rmseS)
corr_all
nanmean(corrS)

figure(2); close(2); figure(2);
set(gcf,'position',[2637,473,1794,829]);
s = 1;
plot(outdata.Mopen(:,s),'linewidth',2); hold on;
%plot(outdata.Moens(:,s),'linewidth',2); hold on;
plot(outdata.Mreg(:,s),'linewidth',2); hold on;
plot(outdata.Val(:,s),'linewidth',2); hold on;
plot(outdata.Obs(:,s),'or');
%plot(outdata.CDF(:,s),'ob');





