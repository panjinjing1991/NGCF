%Code structure
% Copyright(c) Li Lu, 2019

function term = get_terms
    term.get_all_terms = @get_all_terms;
    term.get_depend_terms = @get_depend_terms;
    term.get_lagged_terms = @get_lagged_terms;
    term.get_period_terms = @get_period_terms;
    term.get_spatial_terms = @get_spatial_terms;
    term.get_season_anomaly = @get_season_anomaly;
end

function [depend_terms,lagged_terms,period_terms,spatial_terms,varargout] = ...
          get_all_terms(data_all,...
                        startDate,...
                        lat,lon,...
                        Slen,...
                        nAnnual,nSeasonal,day_lag)
%                   
% exclude leap day of leap year
% 
% Parameters:
% __________
% - data_all: shape as [time,n_feature]
% - startDate: start date of the data
% - lat/lon:
% - Slen:
% - nAnnual/nSeasonal
% - day_lag
%
% Attributes:
% __________
% - depend_terms:
% - lagged_terms:
% - period_terms:
% - spatial_terms:
% - annual_terms:
% - season_terms:
% - unleap_day:
% - jd:
    
% set data
[data_pixel,unleap_day,jd] = remove_leap_day(data_all,startDate,lat,lon);
% get dependent terms
depend_terms = get_depend_terms(data_pixel,day_lag);
% get lagged terms
lagged_terms = get_lagged_terms(data_pixel,day_lag);
% get period terms
[period_terms,annual_terms,season_terms] = get_period_terms...
                                           (data_pixel,...
                                           jd,...
                                           nAnnual,nSeasonal,day_lag);
% get spatial terms
spatial_terms = get_spatial_terms(data_all,Slen,lat,lon,jd,day_lag);
%optional additional outputs
if nargout>4; varargout{1} = annual_terms; end
if nargout>5; varargout{2} = season_terms; end
if nargout>6; varargout{3} = unleap_day; end
if nargout>7; varargout{4} = jd; end

end

function depend_terms = get_depend_terms(data,day_lag)
%
% get dependent terms of input array
% ***attention: only use for our dataframe.(see main.m)

depend_terms = data(day_lag+1:end);
    
end

function [data_,unleap_day,jd] = ...
          remove_leap_day(data,startDate,lat,lon)
%
% exclude leap day of leap year
% 
% Parameters:
% __________
% - data: shape as [time,n_feature]
% - startDate: start date of the data
% - lat/lon:
%
% Attributes:
% __________
% - data_: data with excluding leap day.
% - unleap_day: index of unleap day of year.
% - jd: index of day of year with excluding leap day.
% 
% ***attention: only for Date type like 'yyyy-mm-dd'

% specific data
data = data(:,lon,lat);
% Convert date and time to serial date number
startDateNum = datenum(startDate,'yyyy-mm-dd'); 
% serial date number array
v = (startDateNum:(startDateNum+length(data)-1))'; 
% Convert date and time to vector of components
d = datevec(v); 
%day of year corresponding to data
jd = v - datenum(d(:,1), 1,0);
% day of un-leap year
unleap_day = find(jd~=366); 
N = length(unleap_day);
jd = jd(unleap_day);
data_ = data(unleap_day);

end

function lagged_terms = ...
         get_lagged_terms(data,day_lag)
%
% use extreme index, spatial homeheterogeneity, surface pressure etc.
% to construct X(t-1),X(t-2),...,X(t-dayLag);
%
% Parameters:
% __________
% - data: shape as [time,n_feature]
% - day_lag: lagged day length
%
% Attributes:
% __________
% - lagged_terms: shape as [time,day_lag]

% size of input data
[~,y] = size(data);
% loop for 1:day_lag days
for i = 0:day_lag-1
    % if input data only have one feature, i.e.,soil moisture
    if y==1
        lagged_terms(:,y*i+1:y*(i+1)) = data(i+1:end-(day_lag-i));
    % if input data have more than one feature, i.e.,[soil moisture,press,...]
    else
        lagged_terms(:,y*i+1:y*(i+1)) = data(i+1:end-(day_lag-i),:);
    end
end

end

function [period_terms,annualTerm,seasTerm] = ...
         get_period_terms(data,...
                          jd,...
                          nAnnual,nSeasonal,day_lag)
%     
% Construct periodic terms using Fourier series.
% use method from Tuttle and Savicci,2016,Science.
% 
% Parameters:
% __________
% - nAnnual: 
% - nSeasonal:
% - jd:
% - data:
%
% Attributes:
% __________
% - period_terms:
% - annualTerm:
% - seasTerm:

% size
[N,~] = size(data);
% construct periodic terms varies on annual cycle(18-1.8 year).
annualTerm = zeros(N,nAnnual*2);
d = (1:1:N)';
%for first yrmod term(s), use half sinusoidal cycles 
%(for trend over the 9-yr study period)
for ii = 1
    annualTerm(:,1+(ii-1)*2:ii*2) = [sin(d*(ii)*pi./N),cos(d*(ii)*pi./N)];
end
%for other yrmod terms, use full sinusoidal cycles, 
%starting at 1 cycle per study period
for ii = 1:nAnnual-1
    annualTerm(:,1+((ii+1)-1)*2:(ii+1)*2) = ...
        [sin(d*2*(ii)*pi./N),cos(d*2*(ii)*pi./N)];
end
 
% construct periodic terms varies on seasonal cycle(1year-2.4month).
seasTerm= zeros(N,nSeasonal*2);
for jj = 1:nSeasonal
    seasTerm(:,(jj-1)*2+1:jj*2)=[sin(jd*2*jj*pi./365),cos(jd*2*jj*pi./365)];
end

% % construct lagged P terms using method from Tuttle,2016.
% pLagTerm = zeros(N,2^nPlag-1);
% indexPlag = cell(nPlag+1,1);
% % lag 0 days
% indexPlag(1) = {1:2*(nAnnual+nSeasonal)};
% %need to do 4 times, once for each different # of P lags
% for ii = 1:nPlag
%     
%     pLagLen=2^ii-1; % number of lagged P terms by each day (i.e.,1-4 days)
%     % psub is index of lagged P array
%     %(e.g., 1 for lagged 1 day, 3 for lagged 2 days)
%     if ii == 1
%         psub = 1:pLagLen;
%     elseif ii > 1
%         psub = max(psub)+1:(pLagLen+max(psub));
%     end
%     % index for select independent variable matrix
%     indexPlag(ii+1) = {[1:2*(nAnnual+nSeasonal),2*(nAnnual+nSeasonal)+ psub]};   
%     
%     lagmod = zeros(N,ii);    
%     for jj = 1:ii
%         lagmod(:,jj) = (circshift(POCC,jj));
%     end 
%     
%     % Convert binary number to decimal number, the same decimal number means 
%     % the same condition, give 1 in that particular day, 
%     % +1 for condition that norain in all days
%     JJ = bin2dec(num2str(lagmod))+1; 
%     
%     newlagmod = zeros(N,2^ii);    
%     for mm = 1:N
%         newlagmod(mm,JJ(mm))=1;
%     end
%     
%     % the first column indicate the condition that no rain in all days. 
%     % exclude it
%     pLagTerm(:,psub) = newlagmod(:,2:end);    
% end

% combine period terms and lagged P terms.
period_terms = [annualTerm,seasTerm];
% corresponding to the length of lagged_terms.
period_terms(end-day_lag+1:end,:) = [];
annualTerm(end-day_lag+1:end,:) = [];
seasTerm(1:day_lag,:) = [];
end


function spatial_terms = ...
         get_spatial_terms(data,Slen,lat,lon,jd,day_lag)
%
% get parameters of around pixels, named square. 
% square must be single number
%
% Parameters: 
% __________
% - data: data size as [time,lon,lat] 
%       *** different data mentioned before.
% - Slen: spatial size of the pixels. such as 3x3, i.e., Slen=3
% - lat/lon: selected lat and lon of target pixel
%
% Attributes:
% __________
% - spatial_terms: spatial terms of target pixel.
%                  contains original parameters, maximum-minimun, variance.
%
% ***attention: could add other spatial terms, such as maximum.

if mod(Slen,2)==1
    % get size of data, as [time,lon,lat]
    [Ntime, Nlon, Nlat] = size(data); 
    % initial spatial terms.
    spatial_terms = nan(Ntime,Slen*Slen+1);
    % half of selected square-pixels.
    Hlen = (Slen-1)/2;
    % spatial square range
    lat_range = lat-Hlen:lat+Hlen; 
    lon_range = lon-Hlen:lon+Hlen;
    % numbers of pixels in square 
    N = numel(lat_range)*numel(lon_range);
    % if lat,lon is out of the range
    if sum(lat_range<1)+sum(lon_range<1)==0
        % original parameter of pixel around
        data_ = data(:,lon_range,lat_range);
        data_around = reshape(data_,[Ntime,N]);
        data_pixel = data_around(:,(N-1)/2);
        data_around(:,(N-1)/2) = [];
        data_all = [data_pixel,data_around];
        % construct spatial terms
        spatial_terms = [...
                         data_around,...
                         max(data_all,[],2)-min(data_all,[],2),...
                         sum((data_around-data_pixel).^2,2)...
                         ];
        %              
        spatial_terms = spatial_terms(jd,:);
        spatial_terms(end-day_lag+1:end,:) = [];
    else
        spatial_terms = [];
    end
else
    'square of spatial terms must be odd number'
end
   
end

function extreme_terms = ...
         get_exterme_terms(data,lat,lon)
     
'Need improved'
     
data_pixel = data(:,lon,lat);
data_index = quantile(data_pixel,[0.01,0.05,0.95,0.99]);

% JJ = data_pixel(ii*24+1:ii*24+24);                 
% if length(find(JJ<0.1))>1 & length(find(JJ>0.1))>1
%     extremeIndex(ii,:) = [max(JJ),min(JJ),max(JJ)-min(JJ),std(JJ)...
%         ,length(find(JJ>AA(3))),length(find(JJ>AA(4)))...
%         max(diff([1,find(JJ<0.1)',24]))-1,max([1,find(JJ>0.1)',24])-1];
% elseif length(find(JJ<0.1))==0
%     extremeIndex(ii,:,i,j) = [max(JJ),min(JJ),max(JJ)-min(JJ),std(JJ)...
%         ,length(find(JJ>AA(3))),length(find(JJ>AA(4)))...
%         24,0];
% elseif length(find(JJ>0.1))==0
%     extremeIndex(ii,:,i,j) = [max(JJ),min(JJ),max(JJ)-min(JJ),std(JJ)...
%         ,length(find(JJ>AA(3))),length(find(JJ>AA(4)))...
%         0,24];
% elseif length(find(JJ<0.1))==1
%     extremeIndex(ii,:,i,j) = [max(JJ),min(JJ),max(JJ)-min(JJ),std(JJ)...
%         ,length(find(JJ>AA(3))),length(find(JJ>AA(4)))...
%         max([find(JJ<0.1)-1,24-find(JJ<0.1)]),1];
% elseif length(find(JJ>0.1))==1
%     extremeIndex(ii,:,i,j) = [max(JJ),min(JJ),max(JJ)-min(JJ),std(JJ)...
%         ,length(find(JJ>AA(3))),length(find(JJ>AA(4)))...
%         1,max([find(JJ>0.1)-1,24-find(JJ>0.1)])];
% end

end

function season_anomaly = get_season_anomaly(nSeasonal,...
                                             X,...
                                             season_terms)
%
% seasonal anomaly by method from Tuttle and Savincci[1]
%
% Parameters: 
% __________
% - nSeasonal: 
% - X: target input array
% - season_terms: seasonal terms constructed by get_period_terms()
%
% Attributes:
% __________
% - season_anomaly: input array removed seasonal average,
%

% 
aic_penalty = 2;
N = numel(X);
% index for all possible regression
max_dec_index = 2^nSeasonal-1;
index_ = dec2bin(1:max_dec_index);
index_ = index_ == '1';
% map index, each time remove sinx and cosx.
map_index = nan(2*nSeasonal,1);
for j = 1:nSeasonal
    map_index((j-1)*2+1:(j-1)*2+2) = j;
end
%
const = ones(N,1);
%
aic = nan(max_dec_index,1);
for i = 1:max_dec_index
    II = find(index_(i,:)==1);
    
    JJ = zeros(length(map_index),1);
    for j = 1:length(II)
        JK = find(map_index==II(j));
        JJ(JK) = 1;
    end
    JJ = find(JJ==1);   
    [~,dev] = glmfit([const,season_terms(:,JJ)],X,...
                     'normal',...
                     'link','identity',...
                     'estdisp','on',...
                     'constant','off');
    indicator = index();
    aic(i) = indicator.AIC(numel(JJ)+1,dev,aic_penalty);         
end
%find min aic
[~,ind]=nanmin(aic);
% 
II = find(index_(ind,:)==1);
KK = zeros(length(map_index),1);
for i = 1:length(II)
    JK = find(map_index==II(i));
    KK(JK) = 1;
end
KK = find(KK==1);
% fit best model
[b,~,stats] = glmfit([const,season_terms(:,KK)],X,...
                     'normal',...
                     'link','identity',...
                     'estdisp','on',...
                     'constant','off');
% climatology as determined by best AIC seasonal model of S
clim_ = [const,season_terms(:,KK)]*b;
% subtract mean S climatology from S
season_anomaly = X - clim_;

end
