clear all;clc;
%%
load("Data/inputDATA.mat")

%[p,spImpact,poccPredict,sPredict,Sclim_aic] = spImpactRankOfrandomForestNew(P(:,130,50),S(:,130,63),press(:,130,50));

%%
lon = 130;
lat = 50;
startDate = '2002-06-19';
nAnnual = 6; nSeasonal = 5; day_lag = 4; Slen = 3;
pbcrit=0.1; %threshold to determine light rain and other rain (units: cm)

press = press(1:end-day_lag,:,:);
%%
terms = get_terms();
models = models();
index = index();
%%
[depend_P_terms,lagged_P_terms,annualTerm,seasTerm,spatial_P_terms] = ...
terms.get_all_terms(P,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
[depend_press_terms,lagged_press_terms,~,~,spatial_press_terms] = ...
terms.get_all_terms(press,startDate,lat,lon,Slen,nAnnual,nSeasonal,day_lag);
independ_terms = [lagged_P_terms,lagged_press_terms,...
                  annualTerm,seasTerm];%,...
%                  spatial_P_terms,spatial_press_terms];

POCC = 0.*depend_P_terms;
POCC(find(depend_P_terms>pbcrit))=1;
S = terms.get_all_terms(S,startDate,113-lat,lon,Slen,nAnnual,nSeasonal,day_lag);
%%
[R2P,residP] = models.run_models(independ_terms,POCC,3);
[R2S,residS] = models.run_models(independ_terms,S,3);
%%
season_anomaly = terms.get_season_anomaly(nSeasonal,S,seasTerm);
[impact,p_value] = quatify(S,residS,POCC,residP,season_anomaly);





%%




% No optimized result. adapted RF,ANN,SVR,CNN





% optimized result.




