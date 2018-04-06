clear all
close all
clc

addpath('tools');

N = 1e5;
M = 100;

sig = linspace(0,2.5,M);

B = -50:0.1:50;

for m = 1:M
 T = randn(N,1);
 X = T+randn(N,1)*sig(m);
 [Ixt,Hx,Ht] = mutual_info(X,T,B,B); 
 cc = corrcoef(X,T); rho(m) = cc(2);
 irat(m) = Ixt/Hx;
 snr(m) = std(T)/std(X-T);
 [m/M,sig(m),std(X-T),rho(m),irat(m),snr(m)]
end

figure(1); close(1); fig = figure(1); set(gcf,'color','w')
set(gcf,'position',[10,10,500,300]);

% plot stuff
[ax,h1,h2] = plotyy(sig,[rho;irat],sig,log(snr)); hold on;

set(ax(1),'YColor','k');
set(ax(2),'YColor','k');

h1(1).LineWidth = 3;
h1(2).LineWidth = 3;
h2(1).LineWidth = 3;

h1(1).Color = 'k';
h1(2).Color = 'k';
h2(1).Color = 'k';

h1(1).LineStyle = '-';
h1(2).LineStyle = '--';
h2(1).LineStyle = ':';

grid on;

legend([h1;h2],'linear correlation','information ratio','log(signal/noise)','location','ne');

% labels
xlabel('error standard deviation','fontsize',18);
ylabel(ax(2),'log(signal/noise)','fontsize',18);
ylabel(ax(1),'linear corr. & info rat.','fontsize',18);
title('Error Variance vs. Linear and Nonlinear Statistics','fontsize',16)

fname = 'figures/Figure5_sigVinfo';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



