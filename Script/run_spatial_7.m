clear all; clc; warning('off');
%
load("/work/lilu/NGCF/Data/NCA-LDAS.mat")
load("/work/lilu/NGCF/Data/CONUS_mask.mat")
%
startDate = '2002-06-19';
nAnnual = 6; nSeasonal = 5; day_lag = 4;
Slen = 7; % 3x3:24s, 5x5:29s, 7x7:37s
pbcrit = 0.01;
[nTime,nLon,nLat] = size(P);
% parfor setting
CoreNum=100;
if parpool('local')<=0
    parpool('open','local',CoreNum);
else
    disp('matlab pool already started');
end
% add function attach to parpool
poolobj = gcp;
addAttachedFiles(...
                 poolobj,...
                 {...
                 '/work/lilu/NGCF/NGCF/main.m',...
                 '/work/lilu/NGCF/NGCF/get_terms.m',...
                 '/work/lilu/NGCF/NGCF/models.m',...
                 '/work/lilu/NGCF/NGCF/quatify.m',...
                 '/work/lilu/NGCF/NGCF/index.m',...
                 '/work/lilu/Function/glmfit_itersave.m',...
                 '/work/lilu/Function/mlag.m',...
                 '/work/lilu/Function/statinsertnan.m',...
                 '/work/lilu/Function/statremovenan.m',...
                 '/work/lilu/Function/stattestlink.m',...
                 '/work/lilu/Function/opt_block_length_REV_dec07.m'
                 })
% Question 1.
% The authors used lagged spatial variables (derived from the surrounding 25km 
% pixels) in their model, to account for mesoscale spatial impacts.  However, 
% the variables used were precipitation and pressure.  Air masses can (and do) 
% travel much farther than 75 km in one day, so what good is including these 
% variables in the model? Even if the lag is one day, these variables very 
% likely have no impact on rain occurrence the next day.  If the authors use 
% more than one day lag, then these variables mean even less. So, I have doubts
% that the authors have actually eliminated the impacts of spatial autocorre-
% lation & heterogeneity, as they stated in lines 181-183.
result = nan(5,nLon,nLat);
tic
parfor lon = 1:nLon
    for lat = 1:nLat
        if CONUS_mask_qdeg(lat,lon)==1
            %
            tic
            [impact,p_value,R2P,R2S] = main(P,S,press,...
                                            lon,lat,...
                                            startDate,...
                                            nAnnual,nSeasonal,day_lag,...
                                            Slen,...
                                            pbcrit);
            toc
            result(:,lon,lat) = [impact,p_value,R2P,R2S];
            disp([impact,p_value,R2P,R2S])
        end
    end
end
toc
%
save('/work/lilu/NGCF/Result/result_slen_7.mat','result','-v7.3')
