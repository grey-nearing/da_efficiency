clear all; close all; clc
restoredefaultpath; addpath(genpath(pwd))
%% --- Experiment Parameters ----------------------------------------------

% minimum # of data points to use from each site
Nmin = 1200;

% bins for mutual information
%Nbins = 6;
%bins = logspace(-2.5,-1,Nbins);
bins = [0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1];
Nbins = length(bins);

% number of boostrap samples to use
Nboots = 10;

% total nubmer of samples to test
Nsamps = 20;

% obs/model inexes
oo = 11;
mm = 8;

% init storage
Itot = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Imod = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Iobs = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Icon = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Ifil = zeros(Nbins,Nsamps,Nboots)./0; % info from filter
Efil = zeros(Nbins,Nsamps,Nboots)./0; % efiiciency of filter

%% --- Load Raw Data ------------------------------------------------------

% load data
load('./data/site_data_1.mat');

% dimensions
Is = find(outdata.Ndata(:,oo,mm) >= Nmin);
Ndata = outdata.Ndata(Is,oo,mm);
Nsites = length(Is);
Nall = sum(Nsites);

% init storage
val = zeros(Nall,1); 
mod = zeros(Nall,1); 
obs = zeros(Nall,1); 
fil = zeros(Nall,1); 

% loop through sites
edex = 0;
for s = 1:Nsites; tic

    % load data from sites
    load(strcat('./data/site_data_',num2str(Is(s)),'.mat'));
            
    % pull sites with enough data
    sdex = edex+1; edex = sdex + Ndata(s) - 1;
    val(sdex:edex) = outdata.val(1:Ndata(s));
    mod(sdex:edex) = outdata.Mopen(1:Ndata(s));
    obs(sdex:edex) = outdata.obs(1:Ndata(s));
    fil(sdex:edex) = outdata.Menkf(1:Ndata(s),oo,mm);
    
    % make sure there are no missing values
    assert(min(val)>0); assert(max(val)<=1); assert(~any(isnan(val)));
    assert(min(mod)>0); assert(max(mod)<=1); assert(~any(isnan(mod)));
    assert(min(obs)>0); assert(max(obs)<=1); assert(~any(isnan(obs)));
    assert(min(fil)>0); assert(max(fil)<=1); assert(~any(isnan(fil)));
    
    % screen report
    fprintf('loaded data from site %d of %d; time = %f \n',s,Nsites,toc);
    
end

% total number of data points
Nall = length(val);
assert(Nall==length(mod));
assert(Nall==length(obs));
assert(Nall==length(fil));
fprintf('Total nubmer of data points for this retreival: %d \n',Nall);

% samples
samples = round(logspace(3,log10(Nall),Nsamps));

% make sure there are no missing values
assert(min(val)>0); assert(max(val)<=1); assert(~any(isnan(val)));
assert(min(mod)>0); assert(max(mod)<=1); assert(~any(isnan(mod)));
assert(min(obs)>0); assert(max(obs)<=1); assert(~any(isnan(obs)));
assert(min(fil)>0); assert(max(fil)<=1); assert(~any(isnan(fil)));

%% --- Calculations -------------------------------------------------------

for n = 2:Nbins
    for s = 1:Nsamps
        for b = 1:Nboots
            
            % bins
            B = -0.1:bins(n):1.1;
            
            % pull a sample
            I = datasample(1:Nall,min(samples(s),Nall),'Replace',true);
            X = mod(I); Y = obs(I); Z = val(I); F = fil(I);
            
            %   % only calculate if this is a subsample
            %   if length(I)==Nall && b>1
            %    Imod(n,s,b) = Imod(n,s,1);
            %    Iobs(n,s,b) = Iobs(n,s,1);
            %    Icon(n,s,b) = Icon(n,s,1);
            %    Itot(n,s,b) = Itot(n,s,1);
            %    Ifil(n,s,b) = Ifil(n,s,1);
            %    Efil(n,s,b) = Efil(n,s,1);
            %    continue
            %   end
            
            % screen report
            fprintf('Calculating: bins: %d, samples: %d, bootstrap: %d/%d ... ',length(B),samples(s),b,Nboots); tic;
            
            % total info stats (raw obs)
            [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,Y,Z,B,B,B); % in data
            Imod(n,s,b) = Ixz/Hz;
            Iobs(n,s,b) = Iyz/Hz;
            Icon(n,s,b) = (Iyz-Ixyz)/Hz;
            Itot(n,s,b) = (Ixz+Iyz-Ixyz)/Hz;                        % total available
            
            % total info stats (enkf)
            Ifil(n,s,b) = info(F,Z,B,B,2);                          % info in EnKS
            Efil(n,s,b) = Ifil(n,s,b)/Itot(n,s,b);                  % efficiency of enkf
            Econ(n,s,b) = (Ifil(n,s,b)-Imod(n,s,b))/Icon(n,s,b);    % efficiency of enkf
            
            % screen report
            t = toc; fprintf('finished - time: %d. \n',t);
            
        end % bootstrap
    end % samples
end % bins


%% --- Save Results -------------------------------------------------------

save('./data/2D_convergence_scan.mat','-v7.3');

%% --- Make Figures -------------------------------------------------------

% --- Convergence Plots -------------------------------

for b = 1:Nbins
    binNames(b) = strcat({'bin width = '},num2str(round(bins(b)*1000)/1000));
end

% figure
figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1000,650,1200,1000]);

% plot
subplot(2,1,1);
%errorbar(repmat(samples,Nbins,1)',mean(Itot,3)',-2*std(Itot,[],3)',2*std(Itot,[],3)','-o','linewidth',3);
semilogx(samples,mean(Itot(2:7,:,:),3)','-o','linewidth',2,'markersize',8); hold on;
semilogx(samples,mean(Itot(8:end,:,:),3)','-*','linewidth',2,'markersize',8);
set(gca,'fontsize',18)
grid on;
title('Total Information','fontsize',22)
ylabel('I(Z;X,Y)/H(Z) [nats/nats]','fontsize',20)
xlabel('sample size','fontsize',20);
legend(binNames{2:end},'location',[0.82,0.78,0.1,0.1]);
axis([samples(1),Nall,0,1]);

subplot(2,1,2);
semilogx(samples,mean(Efil(2:7,:,:),3)','-o','linewidth',2,'markersize',8); hold on;
semilogx(samples,mean(Efil(8:end,:,:),3)','-*','linewidth',2,'markersize',8);
set(gca,'fontsize',18)
grid on;
title('Filter Efficiency','fontsize',22)
ylabel('I(Z;X^+)/I(Z;X,Y) [nats/nats]','fontsize',20)
xlabel('sample size','fontsize',20);
% legend(binNames{:},'location','se');
axis([samples(1),Nall,0,1.4]);

fname = strcat('figures/Convergence_MainPlots');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

%% --- Convergence SubPlots -----------------------------

% plot
figure(2); close(2); figure(2);
set(gcf,'color','w','position',[1673        1033        1207         430]);
%errorbar(repmat(samples,Nbins,1)',mean(Itot,3)',-2*std(Itot,[],3)',2*std(Itot,[],3)','-o','linewidth',3);
plot(samples,mean(Itot(2:7,:,:),3)','-o','linewidth',3,'markersize',8); hold on;
plot(samples,mean(Itot(8:end,:,:),3)','-*','linewidth',3,'markersize',8);
set(gca,'fontsize',18)
grid on;
% title('Total Information','fontsize',22)
% ylabel('I(Z;X,Y)/H(Z) [nats/nats]','fontsize',20)
% xlabel('sample size','fontsize',20);
% legend(binNames{:});
axis([samples(1),Nall,0,1]);

fname = strcat('figures/Convergence_SubPlots_Itot');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);

figure(3); close(3); figure(3);
set(gcf,'color','w','position',[1673        1033        1207         430]);
plot(samples,mean(Efil(2:7,:,:),3)','-o','linewidth',3); hold on;
plot(samples,mean(Efil(8:end,:,:),3)','-*','linewidth',3,'markersize',8);
set(gca,'fontsize',18)
grid on;
% title('Filter Efficiency','fontsize',22)
% ylabel('I(Z;X^+)/I(Z;X,Y) [nats/nats]','fontsize',20)
% xlabel('sample size','fontsize',20);
% legend(binNames{:},'location','se');
axis([samples(1),Nall,0,1]);

fname = strcat('figures/Convergence_SubPlots_Efil');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);



