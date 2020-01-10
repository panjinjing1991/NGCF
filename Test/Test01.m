clear all;clc;
%% select nonlinear regression
% contrast different models used default hyperparameters.

% parameters
load('Data/NCA-LDAS.mat')
startDate = '2002-06-23';
nAnnual = 6; nSeasonal = 5; day_lag = 4;
Slen = 3; % 3x3:24s, 5x5:29s, 7x7:37s
pbcrit = 0.01; % thresold for precip
[nTime,nLon,nLat] = size(P); % size for input variable
lon = 100;
lat = 50;
% 0 for only period; 1 plus lag p,press; 
% 2 plus spatial p,press; 3 plus spatial s.
independ_type = 0; 
% 3 for non bias-correct method, 4 for bias-correct method
fun = 4; 

% random 50 point in US
Num = 50;
lon_rand = randperm(nlon,Num);
lat_rand = randperm(nlat,Num);

mains = main();

% initialize
R2_ANN = nan(Num,1);
R2_SVR = nan(Num,1);
R2_RF  = nan(Num,1);

% loop
for i = 1:Num
    i
    lon = lon_rand(i); lat = lat_rand(i);
    [POCC,~,independ_terms,~,~,~,~,P1] = mains.get_sp_terms...
         (P,S,press,lon,lat,startDate,nAnnual,nSeasonal,day_lag,...
          Slen,pbcrit,independ_type);
    if sum(~isnan(P1))~=0
        model = models;
        R2_ANN(i) = model.ANN(independ_terms,POCC);
        R2_SVR(i) = model.SVR(independ_terms,POCC);
        R2_RF(i) = model.RF(independ_terms,P1);
        disp([R2_ANN(i),R2_SVR(i),R2_RF(i)])
    end
    
end

% 
II = find(isnan(R2_ANN));
R2_ANN(II) = [];
R2_RF(II) = [];
R2_SVR(II) = [];
% figure
plot(R2_ANN,'LineWidth',2)
hold on
plot(R2_SVR,'LineWidth',2)
plot(R2_RF,'LineWidth',2)
legend('ANN','SVR','RF')
xlabel('pixel')
ylabel('R2')
axis([1 30 0 1])

