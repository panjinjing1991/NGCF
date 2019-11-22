
load("Data/inputDATA.mat")

lon = 120;
lat = 60;
startDate = '2002-06-19';
nAnnual = 6; nSeasonal = 5; day_lag = 4; Slen = 3;
pbcrit=0.01; %threshold to determine light rain and other rain (units: cm)

% construct terms
[depend_P_terms,lagged_P_terms,period_terms,spatial_P_terms] = ...
          get_terms(P,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
[depend_press_terms,lagged_press_terms,~,spatial_press_terms] = ...
          get_terms(press,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
S = get_terms(S,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);

POCC = 0.*depend_P_terms;
POCC(find(depend_P_terms>pbcrit))=1;

[R2,resid] = models(lagged_P_terms,POCC,3,0.5,0.5);



% No optimized result. adapted RF,ANN,SVR,CNN





% optimized result.




