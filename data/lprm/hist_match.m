function Zm = hist_match(Z,X,dx)

 db = 0.01;
 bins = (0-db/2):db:(1-db/2);

 Hx = hist(X(:),bins);
 Cx = [cumsum(Hx)/sum(Hx)];

 Hz = hist(Z(:),bins);
 Cz = [cumsum(Hz)/sum(Hz)];
 
 Zcdf = interp1(bins,Cz,Z); % Use cdf distribution of Z (i.e., interp1 arguments 1 & 2) to create Z-percentile time series based on the original Z-value time series (i.e., interp1 argument 3)

 try
  Zm = zeros(length(X),1);
  for n = 1:length(X)
   i = find(Zcdf(n)>Cx,1,'last');
   v = bins(i) + db * (Zcdf(n)-Cx(i))/(Cx(i+1)-Cx(i));
   %i = find(Zcdf(n)<=Cx,1,'first');
   %v = bins(i) - db * (Cx(i)-Zcdf(n))/(Cx(i)-Cx(i-1));
   Zm(n) = v;
  end
 catch
  keyboard
 end

%figure(1); hist(Z,bins);
%figure(2); hist(X,bins);
%figure(3); hist(Zm,bins);
%Hm = hist(Zm(:),bins);
%Cm = [cumsum(Hm)/sum(Hm)];
%figure(4);
%plot(Cm,'g')
%hold on
%plot(Cx,'r')
%plot(Cz,'b')
%keyboard

% I1 = find(diff(Xcdf)>0,1,'first');
% I2 = find(diff(Xcdf)>0,1,'last');
% I = I1:I2+1;
%
% try
%  Zm = interp1(Xcdf(I),bins(I),Zcdf_timeseries); % XvalueNew_timeseries is actually "new X-value-similar" or modified/CDF-adjusted Z-value time series 
% catch
%  keyboard
% end
% 


