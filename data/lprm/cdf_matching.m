clear all
close all
clc

Ns = 129;

for s = 1:Ns; tic
    
    % screen report
    fprintf('Running site %d/%d ...',s,Ns);
    
    try
        obs = load(strcat('site_data/lprm_',num2str(s),'.txt'));
        mod = load(strcat('../../open/site_data/open_',num2str(s),'.out'));
        
        assert(size(obs,1) == size(mod,1));
        
        I = find(obs(:,5) > 0);% & obs(:,6) == 6);
        Z = obs(I,5);
        X = mod(I,1);
        Ro = obs(I,6)/100;
        
        if length(I) < 300; error(' '); end
        
        Zc = hist_match(Z,X,0.01);
        In = find(isnan(Zc));
        I(In) = [];
        Zc(In) = [];
        
        Rc = max(0.01, mean(Ro*std(Zc)/std(Z)));
        
        fname = strcat('./site_data/cdf_sig_',num2str(s),'.txt');
        save(fname,'Rc','-ascii');
        
        cdf = obs;
        cdf(:,5) = -9999;
        cdf(I,5) = Zc;
        cdf(:,6) = -9999;
        cdf(I,6) = Rc;
        
        fname = strcat('./site_data/lprm_cdf_',num2str(s),'.txt');
        save(fname,'cdf','-ascii');
        
        fprintf('. finished; time = %f; Ro = %f, Vo = %f Vc = %f, Rc = %f\n',toc,mean(Ro),std(Z),std(Zc),mean(Rc));

    catch
        fprintf('. failed; time = %f \n',toc);
    end
    
end

