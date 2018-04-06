clear all; close all; clc;
restoredefaultpath; addpath(genpath(pwd));
%% --- Experiment Parameters ----------------------------------------------

% number of samples per distribution
N = 1e5;

% number of distributions to test
M = 100;

% root-variances of different distributions to test
sig = linspace(0,2.5,M);

% bins for calculating mutual info
B = -50:0.1:50;

%% --- Do Calcualtions ----------------------------------------------------

for m = 1:M
    T = randn(N,1);
    X = T+randn(N,1)*sig(m);
    [Ixt,Hx,Ht] = mutual_info(X,T,B,B);
    cc = corrcoef(X,T); rho(m) = cc(2);
    irat(m) = Ixt/Hx;
    snr(m) = std(T)/std(X-T);
    [m/M,sig(m),std(X-T),rho(m),irat(m),snr(m)]
end

%% --- Plot Results -------------------------------------------------------

% init figure
figure(1); close(1); fig = figure(1); set(gcf,'color','w')
set(gcf,'position',[10,10,500,300]);

% plot the curves
[ax,h1,h2] = plotyy(sig,[rho;irat],sig,log(snr)); hold on;

% aesthetics
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

% save figure
fname = 'figures/sig_vs_info';
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- END SCRIPT ---------------------------------------------------------


