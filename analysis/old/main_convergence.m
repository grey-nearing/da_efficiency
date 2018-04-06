clear all
close all
clc
restoredefaultpath
addpath(genpath(pwd))

%% *************************************************************************************
%% --- runtime parameters
%% *************************************************************************************

% number of data points to use per site
Nmin = 1200;

% bins for mutual information
%Nbins = 6;
%bins = logspace(-2.5,-1,Nbins);
bins = [0.005,0.01,0.02,0.03,0.04,0.05,0.07];
Nbins = length(bins);

% number of boostrap samples to use
Nboots = 30;

% total nubmer of samples to test
Nsamps = 25;

% init storage
Itot = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Imod = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Iobs = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Icon = zeros(Nbins,Nsamps,Nboots)./0; % total info in obs and model
Ifil = zeros(Nbins,Nsamps,Nboots)./0; % info from filter
Efil = zeros(Nbins,Nsamps,Nboots)./0; % efiiciency of filter

%% *************************************************************************************
%% --- Prep Data
%% *************************************************************************************

% load data
fname = strcat('./data/all_data.mat');
load(fname); 

% dimensions
Ndata = outdata.Ndata;

% separate
val = outdata.Val;
mod = outdata.Moens;
obs = outdata.Obs;
cdf = outdata.CDF;
fil = outdata.Menks;

% pull sites with enough data
I = find(Ndata>=Nmin);
aval = []; amod = []; aobs = []; afil = []; acdf = [];
for s = 1:length(Ndata)
 if Ndata(s) > Nmin
  aval = [aval;val(1:Nmin,s)];
  amod = [amod;mod(1:Nmin,s)];
  aobs = [aobs;obs(1:Nmin,s)];
  afil = [afil;fil(1:Nmin,s)];
  acdf = [acdf;cdf(1:Nmin,s)];
 end
end
val = aval; mod = amod; obs = aobs; fil = afil; cdf = acdf;

% number of sites with data
Nsites = length(find(Ndata>Nmin))

% total number of data points
Nall = length(val)
assert(Nall==length(mod));
assert(Nall==length(obs));
assert(Nall==length(cdf));
assert(Nall==length(fil));
fprintf('Total nubmer of data points for this retreival: %d \n',Nall);

% samples
%samples = round(logspace(3,log10(Nall),Nsamps));
samples = round(linspace(1000,Nall,Nsamps));

% make sure there are no missing values
assert(min(val)>0); assert(max(val)<=1); assert(~any(isnan(val))); 
assert(min(mod)>0); assert(max(mod)<=1); assert(~any(isnan(mod)));
assert(min(obs)>0); assert(max(obs)<=1); assert(~any(isnan(obs)));
assert(min(fil)>0); assert(max(fil)<=1); assert(~any(isnan(fil)));
assert(min(cdf)>0); assert(max(cdf)<=1); assert(~any(isnan(cdf)));

%% *************************************************************************************
%% --- Make SCAN Station Map
%% *************************************************************************************
%make_scan_map()

%% *************************************************************************************
%% --- Calcs
%% *************************************************************************************

for n = 1:Nbins 
 for s = 1:Nsamps
  for b = 1:Nboots

   % bins
   B = -0.1:bins(n):1.1;

   % pull a sample
   I = datasample(1:Nall,min(samples(s),Nall),'Replace',true); 
   X = mod(I); Y = obs(I); Z = val(I); F = fil(I); C = cdf(I);
 
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

   [Ixyz,Ixy,Ixz,Iyz,Hx,Hy,Hz] = mutual_info_3(X,C,Z,B,B,B); % in data
   Cmod(n,s,b) = Ixz/Hz;
   Cobs(n,s,b) = Iyz/Hz;
   Ccon(n,s,b) = (Iyz-Ixyz)/Hz;
   Ctot(n,s,b) = (Ixz+Iyz-Ixyz)/Hz;                        % total available

   % total info stats (enkf)
   Ifil(n,s,b) = info(F,Z,B,B,2);                          % info in EnKS
   Efil(n,s,b) = Ifil(n,s,b)/Itot(n,s,b);                  % efficiency of enkf
   Econ(n,s,b) = (Ifil(n,s,b)-Imod(n,s,b))/Icon(n,s,b);    % efficiency of enkf
   Ecdf(n,s,b) = (Ifil(n,s,b)-Cmod(n,s,b))/Ccon(n,s,b);    % efficiency of enkf

   % screen report
   t = toc; fprintf('finished - time: %d. \n',t);

  end % bootstrap
 end % samples
end % bins


%% *************************************************************************************
%% --- save results 
%% *************************************************************************************
save('./data/2D_convergence.mat');

%% *************************************************************************************
%% --- make plots 
%% *************************************************************************************

mean(Imod(4,end,:))
mean(Iobs(4,end,:))
mean(Icon(4,end,:))
mean(Itot(4,end,:))
1-mean(Ctot(4,end,:)./Itot(4,end,:))
mean(Ifil(4,end,:))
mean(Efil(4,end,:))
mean(Econ(4,end,:))
mean(Ecdf(4,end,:))


% --- Convergence Plots ----------------------------------------------------------------

for b = 1:Nbins
 binNames(b) = strcat({'bin width = '},num2str(round(bins(b)*1000)/1000));
end

figure(1); close(1); figure(1);
set(gcf,'color','w','position',[1000,650,1200,1000]);

% plot
subplot(2,1,1);
errorbar(repmat(samples,Nbins,1)',mean(Itot,3)',-2*std(Itot,[],3)',2*std(Itot,[],3)','-o','linewidth',1); hold on;
plot(samples,mean(Itot,3)','-o','linewidth',2);
grid on;
title('Total Information','fontsize',16)
ylabel('I(Z;X,Y)/H(Z) [nats/nats]','fontsize',16)
xlabel('sample size','fontsize',16);
legend(binNames{:});
axis([1,Nall,0,1]);

subplot(2,1,2);
errorbar(repmat(samples,Nbins,1)',mean(Efil,3)',-2*std(Efil,[],3)',2*std(Efil,[],3)','-o','linewidth',1); hold on;
plot(samples,mean(Efil,3)','-o','linewidth',2);
grid on;
title('Filter Efficiency','fontsize',16)
ylabel('I(Z;X^+)/I(Z;X,Y) [nats/nats]','fontsize',16)
xlabel('sample size','fontsize',16);
legend(binNames{:},'location','se');
axis([1,Nall,0,1]);

fname = strcat('figures/Figure_Convergence');
img = getframe(gcf);
imwrite(img.cdata, [fname, '.png']);


