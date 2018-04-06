clear all
close all
clc

% number of samples
Nsites = 129;
Nobs = 20;
Nmod = 15;

% Data CoubackNdata = zeros(Nsites,Nobs,Nmod)./0;

val = load('../data/scan/site_data/scan_6.txt');
Ntimes = size(val,1);
val_all   = single(zeros(Ntimes,Nsites)./0);
obs_all   = single(zeros(Ntimes,Nsites)./0);
cdf_all   = single(zeros(Ntimes,Nsites)./0);
Mopen_all = single(zeros(Ntimes,Nsites)./0);
Sopen_all = single(zeros(Ntimes,Nsites)./0);

Menkf_all = single(zeros(Ntimes,Nsites,Nobs,Nmod)./0);
Mback_all = single(zeros(Ntimes,Nsites,Nobs,Nmod)./0);
Senkf_all = single(zeros(Ntimes,Nsites,Nobs,Nmod)./0);
Sback_all = single(zeros(Ntimes,Nsites,Nobs,Nmod)./0);

% load sample data
for s = 1:Nsites; ts = tic;

 % load raw data
 try
  val  = single(load(strcat('../data/scan/site_data/scan_'    ,num2str(s),'.txt')));
  obs  = single(load(strcat('../data/lprm/site_data/lprm_'    ,num2str(s),'.txt')));
  cdf  = single(load(strcat('../data/lprm/site_data/lprm_cdf_',num2str(s),'.txt')));
  open = single(load(strcat('../open/site_data/open_'    ,num2str(s),'.out')));
 catch
  fprintf('error loading: Site %d of %d \n',s,Nsites);
  continue
 end

 % extract simultaneous obs
 odex    = find(cdf(:,5)>0);
 vodex   = find(val(odex,5)>0);
 mvodex  = find(open(odex(vodex),1)>0); 

 for o = 1:Nobs
  for m = 1:Nmod; tm = tic;

  % load raw data
   try
    enkf = single(load(strcat('../enkf/site_data/enkf_',num2str(s),'_',num2str(o),'_',num2str(m),'.out')));
   catch
    try
     enkf = single(load(strcat('/discover/nobackup/projects/summa/misc_projects/enkf_storage/enkf_',num2str(s),'_',num2str(o),'_',num2str(m),'.out')));
   catch
     fprintf('error loading enkf: %d, %d, %d ----------------------------- \n',s,o,m);
     continue
    end
   end

   try
    back = single(load(strcat('../enkf/site_data/back_',num2str(s),'_',num2str(o),'_',num2str(m),'.out')));
   catch
    try
     back = single(load(strcat('/discover/nobackup/projects/summa/misc_projects/enkf_storage/back_',num2str(s),'_',num2str(o),'_',num2str(m),'.out')));
   catch
     fprintf('error loading back: %d, %d, %d ----------------------------- \n',s,o,m);
     continue
    end
   end

   % extract simultaneous obs
   emvodex = find(enkf(odex(vodex(mvodex)),1)>0); 
   dex     = odex(vodex(mvodex(emvodex))); 

   % to help throw away sites that don't have enough data
   Ndata(s,o,m) = length(dex);

   if o == 1 && m == 1
    val_all  (1:length(dex),s) = val(dex,5);
    obs_all  (1:length(dex),s) = obs(dex,5);
    cdf_all  (1:length(dex),s) = cdf(dex,5);
    Mopen_all(1:length(dex),s) = open(dex,1);
    Sopen_all(1:length(dex),s) = open(dex,5);
    dexkeep = dex;
   else 
    if o > 1
     assert(Ndata(s,o,m)==Ndata(s,o-1,m))
    end
    if m > 1
     assert(Ndata(s,o,m)==Ndata(s,o,m-1))
    end
    assert(all(dex==dexkeep));
   end
   Menkf_all(1:length(dex),s,o,m) = enkf(dex,1);
   Senkf_all(1:length(dex),s,o,m) = enkf(dex,5);
   Mback_all(1:length(dex),s,o,m) = back(dex,1);
   Sback_all(1:length(dex),s,o,m) = back(dex,5);

   % screen report
   tm = toc(tm);
   fprintf(':        Finished - Site=%d; Obs=%d; Mod=%d; Ndata=%d; time=%f \n',s,o,m,Ndata(s,o,m),tm);

  end %o
 end %s

 % screen report
 ts = toc(ts);
 fprintf('Finished: Site: %d/%d; N = %d; time=%f \n',s,Nsites,Ndata(s,o,m),ts);

end % sites

% ----------------------------------

for s = 1:Nsites

 outdata.Senkf  = squeeze(Senkf_all(:,s,:,:));
 outdata.Sback  = squeeze(Sback_all(:,s,:,:));
 outdata.Menkf  = squeeze(Menkf_all(:,s,:,:));
 outdata.Mback  = squeeze(Mback_all(:,s,:,:));

 outdata.Sopen  = Sopen_all(:,s);
 outdata.Mopen  = Mopen_all(:,s);

 outdata.val    = val_all(:,s);
 outdata.obs    = obs_all(:,s);
 outdata.cdf    = cdf_all(:,s);

 outdata.Ndata  = Ndata;
 outdata.Ndata  = Ndata;

 save(strcat('data/site_data_',num2str(s),'.mat'),'outdata','-v7.3');

end

